function [best_p0, rho] = solve_fixed_time_transfer_indirect(P, t0, tf)
rho = 1;  % Smoothing parameter

% Randomly guess costates until one converges
max_guesses = 100;
iter_tolerance = 1e-5;
for i = 1:max_guesses
    p0_guess = 0.1*rand(7, 1);
    [best_p0, xf_err] = solve_2pbvp(P, t0, tf, rho, p0_guess);
    if norm(xf_err)^2 < iter_tolerance
        break
    end
end

disp('Found initial guess! Continuing with homotopy technique...')

% Use homotopy technique to solve for solution with small rho
for i = 1:5
    rho = rho*0.1;
    [best_p0, xf_err] = solve_2pbvp(P, t0, tf, rho, best_p0);

    if norm(xf_err)^2 > iter_tolerance
        best_p0 = nan;
        return
    end
end