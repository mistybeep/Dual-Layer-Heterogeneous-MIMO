function G = create_G(t, theta, phi, M, N, lambda)
    for m = 1:M
        n_t = [cos(theta(m)) * cos(phi(m)); sin(theta(m))];
        for nt = 1:N
            % 计算信号传播距离差 rho_r
            rho_t = n_t' * t(:, nt);
            % 计算相位差并填入场响应向量
            G(m, nt) = exp(1j * 2 * pi * rho_t / lambda);
        end
    end

end