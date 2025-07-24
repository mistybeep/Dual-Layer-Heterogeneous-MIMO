function W_G = initialize_W_G(H_G, H_L, w_L, gamma, sigma2, K)
    N_G = size(H_G, 2);
    W_G = zeros(N_G, K, 'like', H_G);
    
    % 零强迫初始化
    for k = 1:K
        idx = setdiff(1:K, k);
        H_G_minus_k = H_G(idx, :);
        P_k = eye(N_G) - H_G_minus_k' * pinv(H_G_minus_k * H_G_minus_k') * H_G_minus_k;
        
        h_G_k = H_G(k, :)';
        w_G_k = P_k * h_G_k;
        if norm(w_G_k) < eps
            w_G_k = h_G_k; % 退回到 MRC
        end
        w_G_k = w_G_k / norm(w_G_k);
        
        % 初始功率调整
        signal = abs(H_G(k, :) * w_G_k)^2;
        interference_L = abs(H_L(k, :) * w_L)^2;
        noise = sigma2;
        sinr_required = 2^gamma(k) - 1;
        alpha_k = sqrt(sinr_required * (interference_L + noise) / signal);
        
        W_G(:, k) = alpha_k * w_G_k;
    end
    
    % 考虑所有干扰调整
    % for k = 1:K
    %     signal = abs(H_G(k, :) * W_G(:, k))^2;
    %     interference = sum(abs(H_G(k, :) * W_G(:, setdiff(1:K, k))).^2);
    %     interference_L = abs(H_L(k, :) * w_L)^2;
    %     R_k = log2(1 + signal / (interference + interference_L + sigma2));
    % 
    %     sinr_required = 2^gamma(k) - 1;
    %     sinr_current = signal / (interference + interference_L + sigma2);
    %     alpha = sqrt(sinr_required / max(sinr_current, eps));
    % 
    %     if R_k < gamma(k)
    %         W_G(:, k) = W_G(:, k) * alpha;
    %     elseif R_k > gamma(k) + 0.01
    %         W_G(:, k) = W_G(:, k) / alpha * 0.95; % 平滑缩小
    %     end
    % end
end