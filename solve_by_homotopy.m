function [p0_solution, rho] = solve_by_homotopy(P, t0, tf, rho, p0)
% Solves for a lower-smoothing transfer solution from a higher-smoothing
% solution or other initial guess which is close to optimal

max_failed_iters = 6;
iter_tolerance = 5e-6;
failed_iters = 0;
while failed_iters < max_failed_iters
    [p0_solution, xf_err] = solve_2pbvp(P, t0, tf, rho, p0);

    % Decrease smoothing if the solution converged, otherwise increase it
    if norm(xf_err)^2 > iter_tolerance
        disp('Increasing smoothing...')
        rho = min(rho*10, 1)

        % Perturb intial adjoints slightly
        disturbance = rho/100;
        p0 = p0 + disturbance*(rand(7, 1) - 0.5);

        failed_iters = failed_iters + 1;
    else
        % Done when error and smoothing are below tolerances
        if rho < P.rho_min*10
            disp('Homotopy method complete!')
            return  % Done
        else
            disp('Decreasing smoothing...')
            rho = rho/10
            p0 = p0_solution;
        end
    end
end

disp('Homotopy method failed!')
p0_solution = nan;