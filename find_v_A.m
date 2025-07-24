function v_A = find_v_A(H_GA, H_LA, W_G, w_L, P_J, Sigma_JA, data, sigma2, t_pos_J, r_pos_multi, lambda)

    N_A = size(H_LA, 1);
    theta_JA_t = data.theta_JA_t;
    phi_JA_t = data.phi_JA_t;
    theta_JA_r = data.theta_JA_r;
    phi_JA_r = data.phi_JA_r;
    [H_JA, ~, ~] = generate_channel(t_pos_J, r_pos_multi, lambda, theta_JA_t, phi_JA_t, theta_JA_r, phi_JA_r, Sigma_JA, false);
    v_A = (H_LA*(w_L*w_L')*H_LA' + H_GA*(W_G*W_G')*H_GA' + P_J*(H_JA*H_JA') + sigma2*eye(N_A)) \ (H_LA*w_L);
    v_A = v_A / norm(v_A, 2);

end