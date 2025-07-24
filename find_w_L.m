function w_L = find_w_L(H_G, H_L, H_LA, W_G, v_A, y, z, K, gamma, P, sigma2)

    alpha_k = zeros(K,1);
    for k=1:K
        alpha_k(k) = log2(1 + y(k)) - y(k) ...
                    + 2 * z(k) * sqrt(1 + y(k)) * real(H_G(k,:) * W_G(:, k)) ...
                    - z(k)^2 * (sum_square_abs(H_G(k,:) * W_G) + sigma2) - gamma(k);
    end

    bar_H_LA = H_LA'*(v_A*v_A')*H_LA;
    
    N_L = size(H_L, 2);
    Phi_k = zeros(N_L, N_L, K);
    for k = 1:K
        Phi_k(:, :, k) = (z(k)^2) * (H_L(k, :)' * H_L(k, :));
    end
    Phi = eye(N_L);
    [V, D] = eig(Phi);
    Phi_inv_sqrt = V * diag(1 ./ sqrt(diag(D))) * V'; 

    A = Phi_inv_sqrt * bar_H_LA * Phi_inv_sqrt;

    [eig_vecs, eig_vals] = eig(A);
    eig_vals = diag(eig_vals);
    [~, idx_max] = max(eig_vals);
    u_max = eig_vecs(:, idx_max);
    w_L = sqrt(P) * Phi_inv_sqrt * u_max;
    eta = 0.1*ones(K+1, 1);
    eta2 = 100*ones(K+1, 1);
    iter = 0;
    while(iter < 1000 && norm(eta2 - eta, 2)>1e-4)

        eta2 = eta;
        alpha = P + sum(alpha_k);

        mu = (eta*alpha) / (eta(1)*P + sum(eta(2:K+1) .* alpha_k));
        weighted_Phi2 = zeros(N_L, N_L);
        for k = 1:K
            weighted_Phi2 = weighted_Phi2 + eta(k+1) * Phi_k(:, :, k);
        end
        Phi = mu(1)*eye(N_L) + weighted_Phi2;

        [V, D] = eig(Phi);
        Phi_inv_sqrt = V * diag(1 ./ sqrt(diag(D))) * V';

        A = Phi_inv_sqrt * bar_H_LA * Phi_inv_sqrt;

        [eig_vecs, eig_vals] = eig(A);
        eig_vals = diag(eig_vals);
        [~, idx_max] = max(eig_vals);
        u_max = eig_vecs(:, idx_max);
        w_L = sqrt(alpha) * Phi_inv_sqrt * u_max;
        step_size = 100/(iter+1);
        for i = 1:K+1
            if i == 1
                Phi_temp = eye(N_L);
                alpha_temp = P;
            else
                Phi_temp = Phi_k(:, :, i-1);
                alpha_temp = alpha_k(i-1);
            end
            eta(i) = max(real(eta(i)+step_size*(trace(Phi_temp*(w_L*w_L')) - alpha_temp)),0);
        end
        iter = iter + 1;
    end

end