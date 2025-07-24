% 主脚本

num_simulations = 10; % 仿真次数
rates = zeros(num_simulations, 11); % 存储每次仿真的结果

lambda = 0.01;
norm_A = 4;
M = 8;
N_G = 16;
N_L = 16;
N_J = 16;
N_A = 4;
K = 4;
P_dBm = 0;
P_J_dBm = 30;
sigma2_dBm = -80;
gamma = 1 * ones(K, 1);
err_range = 4 / 180 * pi;
alpha = 2.8;
rho_0_dB = -40;
data = load_dataset('dataset.mat');

% 运行仿真
for sim = 1:num_simulations
    fprintf('运行仿真 %d/%d...\n', sim, num_simulations);
    result = run_simulation(sim, data, norm_A, lambda, M, N_G, N_L, N_J, N_A, K, P_dBm, P_J_dBm, sigma2_dBm, gamma, err_range, alpha, rho_0_dB, Q, R);
    if ~any(isnan(result))
        rates(sim, :) = result;
    else
        rates(sim, :) = [];
    end
end

% 计算平均速率
average_rate = mean(rates, 1);
figure;
plot(0:10,average_rate,'r-o')
ylim([0 max(average_rate)+0.1]);
fprintf('平均速率: %.4f\n', average_rate);

function rate = run_simulation(sim, data, norm_A, lambda, M, N_G, N_L, N_J, N_A, K, P_dBm, P_J_dBm, sigma2_dBm, gamma, err_range, alpha, rho_0_dB, Q, R)
    rng('default')
    P = 10^(P_dBm / 10) * 1e-3; % 转换为 W
    P_J = 10^(P_J_dBm / 10) * 1e-3; % 转换为 W
    sigma2 = 10^(sigma2_dBm / 10) * 1e-3; % 转换为 W

    r_pos_single = zeros(2, K);
    t_G = Gen_UPA(N_G, lambda);
    t_J = Gen_UPA(N_J, lambda);
    t_L = Gen_UPA(N_L, lambda);
    load r_A.mat r_A;
    
    load Sigma.mat Sigma_GA_Test Sigma_GK_Test Sigma_JA_Test Sigma_LA_Test Sigma_LK_Test
    Sigma_GA = squeeze(Sigma_GA_Test(sim,:,:));
    Sigma_LA = squeeze(Sigma_LA_Test(sim,:,:));
    Sigma_JA = squeeze(Sigma_JA_Test(sim,:,:));
    Sigma_GK = squeeze(Sigma_GK_Test(sim, :, :, :));
    Sigma_LK = squeeze(Sigma_LK_Test(sim, :, :, :));
    H_G = zeros(K, N_G);
    for k = 1:K
        r_pos_k = r_pos_single(:, k);
        H_G(k, :) = generate_channel(t_G, r_pos_k, lambda, data.theta_G_t(:, k), data.phi_G_t(:, k), [], [], Sigma_GK(:,:,k), true);
    end
    [H_GA, ~, G_GA] = generate_channel(t_G, r_A, lambda, data.theta_GA_t, data.phi_GA_t, data.theta_GA_r, data.phi_GA_r, Sigma_GA, false);
    H_L = zeros(K, N_L);
    for k = 1:K
        r_pos_k = r_pos_single(:, k);
        H_L(k, :) = generate_channel(t_L, r_pos_k, lambda, data.theta_L_t(:, k), data.phi_L_t(:, k), [], [], Sigma_LK(:,:,k), true);
    end
    [H_LA, ~, G_LA] = generate_channel(t_L, r_A, lambda, data.theta_LA_t, data.phi_LA_t, data.theta_LA_r, data.phi_LA_r, Sigma_LA, false);
    
    [~, ~, G_JA] = generate_channel(t_J, r_A, lambda, data.theta_JA_t, data.phi_JA_t, data.theta_JA_r, data.phi_JA_r, Sigma_JA, false);
    
    w_L = randn(N_L, 1) + 1j * randn(N_L, 1);
    w_L = w_L * sqrt(P) / norm(w_L, 2);
    v_A = randn(N_A, 1) + 1j*randn(N_A, 1);
    v_A = v_A / norm(v_A, 2);
    
    W_G = initialize_W_G(H_G, H_L, w_L, gamma, sigma2, K);
    
    L_max = 10;

    rate_new = cal_SE(r_A, v_A, Sigma_LA, G_LA, Sigma_GA, G_GA, Sigma_JA, G_JA, w_L, W_G, P_J, data, sigma2, lambda);
    rate = rate_new;

    for l = 1:L_max
        y = find_y(H_G, H_L, W_G, w_L, sigma2, K);
    
        z = find_z(H_G, H_L, W_G, w_L, sigma2, y, K);

        v_A = find_v_A(H_GA, H_LA, W_G, w_L, P_J, Sigma_JA, data, sigma2, t_J, r_A, lambda);

        w_L = find_w_L(H_G, H_L, H_LA, W_G, v_A, y, z, K, gamma, P, sigma2);
        
        W_G = find_W_G(H_G, H_GA, H_L, w_L, v_A, y, z, K, gamma, sigma2);
        
        r_A = update_r_A(r_A, data.theta_GA_r, data.phi_GA_r, Sigma_GA, ...
                        data.theta_LA_r, data.phi_LA_r, Sigma_LA, ...
                        data.theta_JA_r, data.phi_JA_r, Sigma_JA, P_J, ...
                        v_A, W_G, w_L, G_GA, G_LA, G_JA, ...
                        lambda, norm_A, sigma2);

        [H_GA, ~, ~] = generate_channel(t_G, r_A, lambda, data.theta_GA_t, data.phi_GA_t, data.theta_GA_r, data.phi_GA_r, Sigma_GA, false);
        [H_LA, ~, ~] = generate_channel(t_L, r_A, lambda, data.theta_LA_t, data.phi_LA_t, data.theta_LA_r, data.phi_LA_r, Sigma_LA, false);

        rate_new = cal_SE(r_A, v_A, Sigma_LA, G_LA, Sigma_GA, G_GA, Sigma_JA, G_JA, w_L, W_G, P_J, data, sigma2, lambda);
        rate = [rate rate_new];
    end
end