function [p0_solution, rho] = solve_fixed_time_transfer_indirect(P, t0, tf)
rho = 1;  % Smoothing parameter

% Randomly guess costates until one converges
max_guesses = 10;
iter_tolerance = 5e-6;
for i = 1:max_guesses
    p0_guess = 0.1*rand(7, 1);
    [best_p0, xf_err] = solve_2pbvp(P, t0, tf, rho, p0_guess);
    if norm(xf_err)^2 < iter_tolerance
        break
    end
end

disp('Found initial guess!')

% Use homotopy technique to solve for solution with small rho
rho = rho*0.1;
p0_solution = solve_by_homotopy(P, t0, tf, rho, best_p0);

end