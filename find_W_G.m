function W_G = find_W_G(H_G, H_GA, H_L, w_L, v_A, y, z, K, gamma, sigma2)
    % 获取维度
    [K, N_G] = size(H_G);
    if isscalar(sigma2)
        sigma2 = sigma2 * ones(K, 1);
    end
    lambda = ones(K,1);

    A = H_GA' * (v_A * v_A') * H_GA;
    for i = 1:K
        h_G_i = H_G(i, :)';
        A = A + lambda(i) * z(i)^2 * (h_G_i * h_G_i');
    end

    W_G_tilde = zeros(N_G, K);
    for k = 1:K
        h_G_k = H_G(k, :)';
        w_G_k = pinv(A) * h_G_k;
        W_G_tilde(:, k) = w_G_k / norm(w_G_k, 2);
    end

    c1 = zeros(K, 1);
    c2 = zeros(K, K);
    beta = zeros(K, 1);
    for k = 1:K
        h_G_k = H_G(k, :)';
        h_L_k = H_L(k, :)';
        c1(k) = abs(h_G_k' * W_G_tilde(:, k))^2;
        for i = 1:K
            c2(k, i) = (2^gamma(k) - 1)*abs(h_G_k' * W_G_tilde(:, i))^2;
        end
        beta(k) = (2^gamma(k) - 1)*(abs(h_L_k' * w_L)^2 + sigma2(k));
    end

    C = diag(c1);
    for k = 1:K
        for i = 1:K
            if i ~= k
                C(k, i) = -c2(k, i);
            end
        end
    end

    p = C \ beta;

    W_G = W_G_tilde * diag(sqrt(p));

end