function [tf_list, dv_list, p0_list] = build_pareto_front(P, t0, init_tf, ...
    rho, init_p0, tf_min, tf_max, dtf)

tf = init_tf;
p0 = init_p0;
tf_list = [];
dv_list = [];
p0_list = [];

% Loop, increasing tf and re-solving each time
while tf <= tf_max
    % Solve the problem and propagate the solution
    p0 = solve_by_homotopy(P, t0, tf, rho, p0);
    
    % If no solution was found continue to the next iteration
    if isnan(p0)
        continue
    end

    [~, X] = propagator_MEE_indirect(P, t0, tf, p0(1:6), p0(7), rho);

    % Calculate delta-V
    dv = calc_delta_v(P, X);

    % Record transfer time, delta-V, and initial adjoints
    tf_list = [tf_list, tf];
    dv_list = [dv_list, dv];
    p0_list = [p0_list, p0];
    tf = tf + dtf
end

% Reset transfer time and adjoints
tf = init_tf;
p0 = init_p0;

% Loop, decreasing tf and re-solving each time
while tf >= tf_min
    % Solve the problem and propagate the solution
    p0 = solve_by_homotopy(P, t0, tf, rho, p0);
    
    % If no solution was found continue to the next iteration
    if isnan(p0)
        continue
    end

    [~, X] = propagator_MEE_indirect(P, t0, tf, p0(1:6), p0(7), rho);

    % Calculate delta-V
    dv = calc_delta_v(P, X);

    % Record transfer time, delta-V, and initial adjoints
    tf_list = [tf_list, tf];
    dv_list = [dv_list, dv];
    p0_list = [p0_list, p0];
    tf = tf - dtf
end