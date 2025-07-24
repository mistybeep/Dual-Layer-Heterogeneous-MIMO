function F = create_F(r_A, theta, phi, M, N, lambda)
    for m = 1:M
        n_r = [cos(theta(m)) * cos(phi(m)); cos(theta(m)) * sin(phi(m))];
        for nr = 1:N
            % 计算信号传播距离差 rho_r
            rho_r = n_r' * r_A(:, nr);
            % 计算相位差并填入场响应向量
            F(m, nr) = exp(1j * 2 * pi * rho_r / lambda);
        end
    end

end