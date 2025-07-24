clear all;
clc;
close all;

rng(42,'twister');
lambda = 0.01;
norm_A = 3;
M = 8;
N_G = 8;
N_L = 8;
N_J = 8; 
N_A = 4;
d = 40;
K = 4;
P_dBm = 20; 
P = 10^(P_dBm / 10) * 1e-3;
P_J = 10^(30 / 10) * 1e-3;
sigma2_dBm = -80;
sigma2 = 10^(sigma2_dBm / 10) * 1e-3; % 转换为线性单位 (W)
gamma = 1*ones(K, 1);
err_range = 4/180 * pi;
alpha = 2.6;
rho_0_dB = -30;
rho_0 = 10^(rho_0_dB / 10);

Q = 21;
R = 21;

user_center = [40; 30; 0]; % 用户区域中心
radius = 10; % 用户区域半径
GBS_pos = [0; 0; 30]; % GBS 位置
LBS_pos = [0; 100; 30]; % LBS 位置
UAV_pos = [10; 60; 100]; % 无人机位置
Jammer_pos = [40; 60; 0]; % 干扰源位置

% t_G = generate_upa_coordinates(norm_A, N_G, lambda);
% t_L = generate_upa_coordinates(norm_A, N_L, lambda);
% t_J = generate_upa_coordinates(norm_A, N_J, lambda);
% r_pos_single = zeros(2, K); % 单天线用户位置（假设在原点）
% r_A = generate_upa_coordinates(norm_A, N_A, lambda);
load('coord.mat')
load('dataset.mat')

% 生成路径响应矩阵 Sigma (对角矩阵)
Sigma_GA = generate_Sigma(M, rho_0, alpha, d_GA);
% Sigma_GA = Sigma_GA*1e3;
Sigma_LA = generate_Sigma(M, rho_0, alpha, d_LA);
Sigma_JA = generate_Sigma(M, rho_0, alpha, d_JA);
% load('Sigma.mat')
% 生成基站 G 的信道
H_G = zeros(K, N_G); % 基站 G 到 K 个单天线用户的信道 (K x Nt)
Sigma_GK = zeros(M, M, K);
G_GK = zeros(M, N_G, K);
for k = 1:K
    r_pos_k = r_pos_single(:, k); % 第 k 个单天线用户的位置 (2 x 1)
    Sigma_GK(:,:,k) = generate_Sigma(M, rho_0, alpha, d_single_G(k));
    [H_G(k, :),~,G_k] = generate_channel(t_G, r_pos_k, lambda, theta_G_t(:, k), phi_G_t(:, k), [], [], Sigma_GK(:,:,k), true);
    G_GK(:,:,k) = G_k;
    % H_G(k, :) = H_G(k, :)*1e3;
    % Sigma_GK(:,:,k) = Sigma_GK(:,:,k)*1e3;
end
[H_GA, F_GA, G_GA] = generate_channel(t_G, r_A, lambda, theta_GA_t, phi_GA_t, theta_GA_r, phi_GA_r, Sigma_GA, false);
% H_GA = H_GA*1e3;
% 生成基站 L 的信道
H_L = zeros(K, N_L); % 基站 L 到 K 个单天线用户的信道 (K x Nt)
Sigma_LK = zeros(M, M, K);
for k = 1:K
    r_pos_k = r_pos_single(:, k); % 第 k 个单天线用户的位置 (2 x 1)
    Sigma_LK(:,:,k) = generate_Sigma(M, rho_0, alpha, d_single_L(k));
    H_L(k, :) = generate_channel(t_L, r_pos_k, lambda, theta_L_t(:, k), phi_L_t(:, k), [], [], Sigma_LK(:,:,k), true);
    % H_L(k, :) = H_L(k, :)*1e3;
    % Sigma_LK(:,:,k) = Sigma_LK(:,:,k)*1e3;
end
[H_LA, F_LA, G_LA] = generate_channel(t_L, r_A, lambda, theta_LA_t, phi_LA_t, theta_LA_r, phi_LA_r, Sigma_LA, false);
% H_LA = H_LA*1e3;
JA = struct();
JA.N = N_J;
JA.M = M;
JA.theta_t = theta_JA_t;
JA.phi_t = phi_JA_t;
JA.theta_r = theta_JA_r;
JA.phi_r = phi_JA_r;
JA.Sigma = Sigma_JA;


% save("Sigma.mat", 'Sigma_GK', 'Sigma_LK', 'Sigma_GA', 'Sigma_LA', ...
%          'Sigma_JA');
% sigma2 = sigma2*1e6;

[H_JA, F_JA, G_JA] = generate_channel(t_J, r_A, lambda, theta_JA_t, phi_JA_t, theta_JA_r, phi_JA_r, Sigma_JA, false);

% W_G = H_G';
% W_G = H_G' * inv(H_G * H_G');
% H_L_null = null(H_L); % 计算 H_L 的零空间
w_L = randn(N_L, 1) + 1j * randn(N_L, 1);
v_A = randn(N_A, 1) + 1j*randn(N_A, 1);
% save("beam.mat", 'w_L', 'v_A')
% load("beam.mat");
w_L = w_L * sqrt(P) / norm(w_L, 2);

v_A = v_A / norm(v_A, 2);

W_G = initialize_W_G(H_G, H_L, w_L, gamma, sigma2, K);

mu = ones(Q,1)/Q;
nu = ones(R,1)/R;
sum_A = zeros(N_A, N_A);
for q = 1:Q
    theta_q = theta_JA_r - err_range/2 + (q-1)*err_range/(Q-1);
    for r = 1:R
        phi_r = phi_JA_r - err_range/2 + (r-1)*err_range/(R-1);
        F_JA = create_F(r_A, theta_q, phi_r, M, N_A, lambda);
        sum_A = sum_A + mu(q)*nu(r)*F_JA'*Sigma_JA*(G_JA*G_JA')*Sigma_JA'*F_JA;
    end
end
kappa_2 = real(abs(v_A'*F_LA'*Sigma_LA*G_LA*w_L)^2 / ...
    real(sum(abs(v_A' * F_GA'* Sigma_GA * G_GA * W_G).^2) ...
    + P_J*v_A'*sum_A*v_A + norm(v_A, 2)^2*sigma2));
rate_new = log2(1+kappa_2);
L_max = 20;

% z = ones(K,1);
% y = ones(K,1);
% 
% [gradient, bar_g_G] = cal_gradient_L1(t_G, theta_GA_t, phi_GA_t, Sigma_GA, G_GA, F_GA, ...
%                                                 theta_G_t, phi_G_t, Sigma_GK, G_GK, ...
%                                                 z, y, v_A, W_G, ...
%                                                 20, lambda, 2);
% [hessian, max_eigenvalue] = cal_hessian_L1(t_G, theta_GA_t, phi_GA_t, Sigma_GA, F_GA, ...
%                                                     theta_G_t, phi_G_t, Sigma_GK, ...
%                                                     z, y, v_A, W_G, bar_g_G, ...
%                                                     20, lambda, 2);
% 
rate = rate_new;

for l = 1:L_max
    y = find_y(H_G, H_L, W_G, w_L, sigma2, K);

    z = find_z(H_G, H_L, W_G, w_L, sigma2, y, K);
    % w_L = find_w_L(H_G, H_L, H_LA, W_G, v_A, y, z, K, gamma, P, sigma2);
    w_L = find_w_L_2(H_G, H_L, H_LA, W_G, v_A, y, z, K, gamma, P, sigma2);

    % W_G = find_W_G(H_G, H_GA, H_L, w_L, v_A, y, z, K, gamma, sigma2);
    % norm(W_G, 'fro')
    % W_G = find_W_G_2(H_G, H_GA, H_L, w_L, v_A, y, z, K, gamma, sigma2);
    W_G = find_W_G_3(H_G, H_GA, H_L, w_L, v_A, y, z, K, gamma, sigma2);

    [v_A, mu, nu, A] = find_v_A(H_GA, H_LA, W_G, w_L, v_A, P_J, Q, R, JA, err_range, sigma2, t_J, r_A, lambda);
    
    sum_A = zeros(N_A, N_A);    
    for q = 1:Q
        for r = 1:R
            sum_A = sum_A + mu(q)*nu(r)*A(:,:,q,r);
        end
    end
    % real(abs(v_A'*H_LA*w_L)^2 / real(sum(abs(v_A' * H_GA * W_G).^2) + P_J*v_A'*sum_A*v_A + norm(v_A, 2)^2*sigma2))
    % r_A = update_r_A(r_A, theta_GA_r, phi_GA_r, Sigma_GA, ...
    %                         theta_LA_r, phi_LA_r, Sigma_LA, ...
    %                         theta_JA_r, phi_JA_r, Sigma_JA, P_J, ...
    %                         v_A, W_G, w_L, G_GA, G_LA, G_JA, ...
    %                         mu, nu, lambda, err_range, norm_A, sigma2);
    % 
    % [H_GA, ~, ~] = generate_channel(t_G, r_A, lambda, theta_GA_t, phi_GA_t, theta_GA_r, phi_GA_r, Sigma_GA, false);
    % [H_LA, ~, ~] = generate_channel(t_L, r_A, lambda, theta_LA_t, phi_LA_t, theta_LA_r, phi_LA_r, Sigma_LA, false);
    % % [H_JA, ~, ~] = generate_channel(t_J, r_A, lambda, theta_JA_t, phi_JA_t, theta_JA_r, phi_JA_r, Sigma_JA, false);
    % F_GA = create_F(r_A, theta_GA_r, phi_GA_r, M, N_A, lambda);
    % F_LA = create_F(r_A, theta_LA_r, phi_LA_r, M, N_A, lambda);
    % sum_A = zeros(N_A, N_A);
    % for q = 1:Q
    %     theta_q = theta_JA_r - err_range/2 + (q-1)*err_range/(Q-1);
    %     for r = 1:R
    %         phi_r = phi_JA_r - err_range/2 + (r-1)*err_range/(R-1);
    %         F_JA = create_F(r_A, theta_q, phi_r, M, N_A, lambda);
    %         sum_A = sum_A + mu(q)*nu(r)*F_JA'*Sigma_JA*(G_JA*G_JA')*Sigma_JA'*F_JA;
    %     end
    % end
    % kappa_2 = real(abs(v_A'*F_LA'*Sigma_LA*G_LA*w_L)^2 / ...
    %     real(sum(abs(v_A' * F_GA'* Sigma_GA * G_GA * W_G).^2) ...
    %     + P_J*v_A'*sum_A*v_A + norm(v_A, 2)^2*sigma2));
    % rate_new = log2(1+kappa_2);

    % t_G = update_t_G(t_G, theta_GA_t, phi_GA_t, theta_G_t, phi_G_t, ...
    %                       Sigma_GA, Sigma_GK, F_GA, v_A, W_G, ...
    %                       H_L, w_L, y, z, gamma, norm_A, lambda, sigma2);
    t_G = update_t_G(t_G, theta_GA_t, phi_GA_t, theta_G_t, phi_G_t, ...
                            Sigma_GA, Sigma_GK, F_GA, v_A, W_G, ...
                            H_L, w_L, y, z, gamma, norm_A, lambda, sigma2);

    for k = 1:K
        r_pos_k = r_pos_single(:, k);
        H_G(k, :) = generate_channel(t_G, r_pos_k, lambda, theta_G_t(:, k), phi_G_t(:, k), [], [], Sigma_GK(:,:,k), true);
        % H_G(k, :) = generate_channel(t_G, r_pos_k, lambda, theta_G_t(:, k), phi_G_t(:, k), [], [], Sigma_GK(:,:,k), true);
    end
    % [H_GA, F_GA, G_GA] = generate_channel(t_G, r_A, lambda, theta_GA_t, phi_GA_t, theta_GA_r, phi_GA_r, Sigma_GA, false);
    [H_GA, ~, ~] = generate_channel(t_G, r_A, lambda, theta_GA_t, phi_GA_t, theta_GA_r, phi_GA_r, Sigma_GA, false);
    % 
    rate_new = log2(1+real(abs(v_A'*H_LA*w_L)^2 / real(sum(abs(v_A' * H_GA * W_G).^2) + P_J*v_A'*sum_A*v_A + norm(v_A, 2)^2*sigma2)));
    rate = [rate rate_new];

end

for k = 1:K
    h_G_k = H_G(k, :)';
    Gamma_k = sum(abs(h_G_k' * W_G).^2) + abs(H_L(k,:) * w_L).^2 + sigma2;
    % r(k) = log2(1+abs(h_G_k' * W_G(:, k))^2 / (Gamma_k - abs(h_G_k' * W_G(:, k))^2));
    r(k) = log2(1+y(k))- y(k) + (1+y(k))*abs(h_G_k' * W_G(:, k))^2/Gamma_k;
end
r
plot(0:L_max,rate,'r-o');
ylim([0 max(rate)+0.1]);