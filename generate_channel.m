function [H, F, G] = generate_channel(t_pos, r_pos, lambda, theta_t, phi_t, theta_r, phi_r, Sigma, is_single_rx)
    % 输入参数:
    %   t_pos: 发射端天线位置矩阵 (2 x Nt)，每列为 [x_t; y_t]
    %   r_pos: 接收端天线位置矩阵 (2 x Nr)，每列为 [x_r; y_r]，单天线时建议为 [0;0]
    %   M: 路径数（发射端和接收端相同）
    %   lambda: 载波波长
    %   theta_t: 发射端仰角 AoD (M x 1)
    %   phi_t: 发射端方位角 AoD (M x 1)
    %   theta_r: 接收端仰角 AoA (M x 1)，单天线时忽略
    %   phi_r: 接收端方位角 AoA (M x 1)，单天线时忽略
    %   is_single_rx: 布尔值，true 表示接收端单天线，false 表示多天线
    %   d: 发射端到接收端的距离 (标量，单位：米)
    % 输出:
    %   H: 信道矩阵 (Nr x Nt)，单天线时 Nr = 1

    Nt = size(t_pos, 2); % 发射端天线数
    Nr = size(r_pos, 2); % 接收端天线数
    M = size(theta_t, 1);

    % 初始化发射端场响应矩阵
    G = zeros(M, Nt); % 发射端场响应矩阵 G (M x Nt)

    % 计算发射端波矢量和场响应向量
    for m = 1:M
        % 发射端归一化波矢量 n_t
        n_t = [cos(theta_t(m)) * cos(phi_t(m)); sin(theta_t(m))];
        for nt = 1:Nt
            % 计算信号传播距离差 rho_t
            rho_t = n_t' * t_pos(:, nt);
            % 计算相位差并填入场响应向量
            G(m, nt) = exp(1j * 2 * pi * rho_t / lambda);
        end
    end

    % 初始化接收端场响应矩阵
    F = zeros(M, Nr); % 接收端场响应矩阵 F (M x Nr)

    % 根据接收端是否为单天线计算 F
    if is_single_rx
        % 单天线情况：F 为全 1 向量 (M x 1)
        if Nr ~= 1
            error('单天线模式下接收端天线数必须为 1');
        end
        F = ones(M, Nr); % Nr = 1 时，F 为 M x 1 全 1 向量
    else
        % 多天线情况：根据角度计算 F
        if nargin < 9 || isempty(theta_r) || isempty(phi_r)
            error('多天线模式下需要提供 theta_r 和 phi_r');
        end
        for m = 1:M
            % 接收端归一化波矢量 n_r
            n_r = [cos(theta_r(m)) * cos(phi_r(m)); cos(theta_r(m)) * sin(phi_r(m))];
            for nr = 1:Nr
                % 计算信号传播距离差 rho_r
                rho_r = n_r' * r_pos(:, nr);
                % 计算相位差并填入场响应向量
                F(m, nr) = exp(1j * 2 * pi * rho_r / lambda);
            end
        end
    end 

    % 计算信道矩阵 H = F^H * Sigma * G
    H = F' * Sigma * G;

end