function p0_solution = solve_by_homotopy(P, t0, tf, rho, p0)

max_iters = 15;
iter_tolerance = 5e-6;
for i = 1:max_iters
    [p0_solution, xf_err] = solve_2pbvp(P, t0, tf, rho, p0);
    
    if rho <= P.rho_min
        disp('Homotopy method complete!')
        return  % Done
    end

    % Decrease smoothing if the converged, otherwise increase it
    if norm(xf_err)^2 > iter_tolerance
        disp('Increasing smoothing...')
        rho = rho*10
    else
        disp('Decreasing smoothing...')
        rho = rho*0.1
        p0 = p0_solution;
    end
end

disp('Homotopy method failed!')
p0_solution = nan;