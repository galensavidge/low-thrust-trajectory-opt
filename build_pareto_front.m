clear
clc
close all
format long g

addpath('orbit_util')
addpath('orbit_util\autogen\')

% Set up the problem
P = problem_setup();

t0 = 0;
day = 86400/P.TU;  % One day in normalized time units
init_tf = 180*day;  % Start with a six month transfer

% Solve one transfer fully to get an initial guess
[init_p0, rho_init] = solve_fixed_time_transfer_indirect(P, t0, init_tf);

% Pareto front time bounds
tf_min = 160*day/P.TU;
tf_max = 220*day/P.TU;
dtf = 2.5*day/P.TU;

tf = init_tf;
p0 = init_p0;
tf_list = tf;
[~, X] = propagator_MEE_indirect(P, t0, tf, p0(1:6), p0(7), rho_init);
dv_list = calc_delta_v(P, X);
p0_list = init_p0;

% Loop, increasing tf and re-solving each time
while tf <= tf_max
    tf = tf + dtf

    % Solve the problem and propagate the solution
    [p0_solution, rho_solution] = solve_by_homotopy(P, t0, tf, rho_init*10, p0);
    
    % If no solution was found continue to the next iteration
    if isnan(p0_solution)
        continue
    end

    p0 = p0_solution;

    [~, X] = propagator_MEE_indirect(P, t0, tf, p0(1:6), p0(7), rho_solution);

    % Calculate delta-V
    dv = calc_delta_v(P, X);

    % Record transfer time, delta-V, and initial adjoints
    tf_list = [tf_list, tf];
    dv_list = [dv_list, dv];
    p0_list = [p0_list, p0];
end

% Reset transfer time and adjoints
tf = init_tf;
p0 = init_p0;

% Loop, decreasing tf and re-solving each time
while tf >= tf_min
    tf = tf - dtf

    % Solve the problem and propagate the solution
    [p0_solution, rho_solution] = solve_by_homotopy(P, t0, tf, rho_init*10, p0);
    
    % If no solution was found continue to the next iteration
    if isnan(p0_solution)
        continue
    end

    p0 = p0_solution;

    [~, X] = propagator_MEE_indirect(P, t0, tf, p0(1:6), p0(7), rho_solution);

    % Calculate delta-V
    dv = calc_delta_v(P, X);

    % Record transfer time, delta-V, and initial adjoints
    tf_list = [tf_list, tf];
    dv_list = [dv_list, dv];
    p0_list = [p0_list, p0];
end

% Save results to a file
filename = sprintf('optimal_transfers_%s', datestr(now,'mm-dd-yyyy HH-MM'));
save(filename, 'tf_list', 'dv_list', 'p0_list')