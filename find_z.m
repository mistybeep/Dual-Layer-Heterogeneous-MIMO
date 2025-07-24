function z = find_z(H_G, H_L, W_G, w_L, sigma2, y, K)
    % numerator = sqrt(1+y).*real(diag(H_G * W_G));
    % denominator = sum(abs(H_G * W_G).^2, 2) + abs(H_L * w_L).^2 + sigma2*ones(K,1);
    % z = numerator./denominator; + abs(H_L(k,:) * w_L)^2
    z = ones(K, 1); % 初始 z
    for k = 1:K
        h_G_k = H_G(k, :)';
        Gamma_k = sum(abs(h_G_k' * W_G).^2) + abs(H_L(k,:) * w_L)^2 + sigma2;
        % y(k) = abs(h_G_k' * W_G(:, k))^2 / (Gamma_k - abs(h_G_k' * W_G(:, k))^2);
        z(k) = sqrt(1 + y(k)) * real(h_G_k' * W_G(:, k)) / Gamma_k;
    end

end