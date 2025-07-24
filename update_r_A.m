function r_A = update_r_A(r_A, theta_GA, phi_GA, Sigma_GA, ...
                            theta_LA, phi_LA, Sigma_LA, ...
                            theta_JA, phi_JA, Sigma_JA, P_J, ...
                            v_A, W_G, w_L, G_GA, G_LA, G_JA, ...
                            lambda, normA, sigma2)

    N_A = size(v_A, 1);
    M_GA = size(theta_GA, 1);
    M_LA = size(theta_LA, 1);
    M_JA = size(theta_JA, 1);

    t = 0;
    max_iter = 100;
    tol = 1e-3;
    kappa = 0;
    kappa2 = -100;

    while(abs(kappa - kappa2)>tol && t<max_iter)
        t = t + 1;

        kappa2 = kappa;
        F_GA = create_F(r_A, theta_GA, phi_GA, M_GA, N_A, lambda);
        F_LA = create_F(r_A, theta_LA, phi_LA, M_LA, N_A, lambda);
        F_JA = create_F(r_A, theta_JA, phi_JA, M_JA, N_A, lambda);
        H_JA = F_JA'*Sigma_JA*G_JA;
        kappa = real(abs(v_A'*F_LA'*Sigma_LA*G_LA*w_L)^2 / ...
                     real(sum(abs(v_A' * F_GA'* Sigma_GA * G_GA * W_G).^2) ...
                     + P_J*v_A'*(H_JA*H_JA')*v_A + norm(v_A, 2)^2*sigma2));

        for n = 1:N_A
            F_GA = create_F(r_A, theta_GA, phi_GA, M_GA, N_A, lambda);
            F_LA = create_F(r_A, theta_LA, phi_LA, M_LA, N_A, lambda);
            

            [gradient, bar_f_A, bar_f_JA] = cal_gradient(r_A, theta_GA, phi_GA, Sigma_GA, ...
                                                            theta_LA, phi_LA, Sigma_LA, F_GA, F_LA, ...
                                                            theta_JA, phi_JA, Sigma_JA, P_J, ...
                                                            v_A, W_G, w_L, G_GA, G_LA, G_JA, kappa, n, lambda);
            [~, tau] = cal_hessian(r_A, theta_GA, phi_GA, Sigma_GA, ...
                                    theta_LA, phi_LA, Sigma_LA, ...
                                    theta_JA, phi_JA, Sigma_JA, P_J, ...
                                    v_A, W_G, w_L, G_GA, G_LA, G_JA, bar_f_A, bar_f_JA, ...
                                    kappa, n, lambda);
            
            r_Astar = r_A(:, n) - real(gradient) / real(tau);
            if n==1
                r_A_opt = r_A(:, 2:end);
                r_An = find_pos(r_A_opt.', r_Astar.', lambda, normA, norm(real(gradient) / real(tau), 2));
                if isnan(r_An(1,1)) || isnan(r_An(1,2))
                    r_A(:, n) = r_A(:, n);
                else
                    r_A(:, n) = r_An.';
                end
            else
                r_A_opt = r_A(:, 1:n-1);
                if n<N_A
                    r_A_opt = [r_A_opt, r_A(:, n+1:end)];
                end
                r_An = find_pos(r_A_opt.', r_Astar.', lambda, normA, norm(real(gradient) / real(tau), 2));
                if isnan(r_An(1,1)) || isnan(r_An(1,2))
                    r_A(:, n) = r_A(:, n);
                else
                    r_A(:, n) = r_An.';
                end
            end

        end
    end

end