function y = find_y(H_G, H_L, W_G, w_L, sigma2, K)
    % numerator = abs(diag(H_G * W_G)).^2;
    % interference = sum(abs(H_G * W_G).^2, 2) - numerator;
    % interference_LBS = abs(H_L * w_L).^2;
    % denominator = interference + interference_LBS + sigma2*ones(K,1); + abs(H_L(k,:) * w_L)^2
    % 
    % y = numerator ./ denominator;
    y = zeros(K, 1); % 初始 y
    
    for k = 1:K
        h_G_k = H_G(k, :)';
        Gamma_k = sum(abs(h_G_k' * W_G).^2) + abs(H_L(k,:) * w_L)^2 + sigma2;
        y(k) = abs(h_G_k' * W_G(:, k))^2 / (Gamma_k - abs(h_G_k' * W_G(:, k))^2);
        % z(k) = sqrt(1 + y(k)) * real(h_G_k' * W_G(:, k)) / Gamma_k;
    end

end