function SE = cal_SE(r_A, v_A, Sigma_LA, G_LA, Sigma_GA, G_GA, Sigma_JA, G_JA, w_L, W_G, P_J, data, sigma2, lambda)
    N_A = size(v_A, 1);
    M_JA = size(Sigma_JA, 1);
    M_GA = size(Sigma_GA, 1);
    M_LA = size(Sigma_LA, 1);

    theta_JA_r = data.theta_JA_r;
    phi_JA_r = data.phi_JA_r;

    F_JA = create_F(r_A, theta_JA_r, phi_JA_r, M_JA, N_A, lambda);  % 调用 create_F 函数
    H_JA = F_JA' * Sigma_JA * G_JA;

    % 计算 F_GA 和 F_LA
    F_GA = create_F(r_A, data.theta_GA_r, data.phi_GA_r, M_GA, N_A, lambda);
    F_LA = create_F(r_A, data.theta_LA_r, data.phi_LA_r, M_LA, N_A, lambda);

    % 计算分子
    numerator = abs(v_A' * F_LA' * Sigma_LA * G_LA * w_L)^2;

    % 计算分母的三项
    term1 = sum(abs(v_A' * F_GA' * Sigma_GA * G_GA * W_G).^2);
    term2 = P_J * (v_A' * (H_JA * H_JA') * v_A);
    term3 = norm(v_A, 2)^2 * sigma2;
    denominator = term1 + term2 + term3;

    % 计算 SINR
    sinr = real(numerator / denominator);  % 取实部

    % 计算频谱效率 (SE)
    SE = log2(1 + sinr);
end