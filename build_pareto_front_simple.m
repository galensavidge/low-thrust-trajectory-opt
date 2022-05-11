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

% Pareto front time bounds
tf_min = 160*day;
tf_max = 220*day;
dtf = 2.5*day;

tf_list = [];
dv_list = [];
p0_list = [];

tf = tf_min

% Loop, increasing tf and re-solving each time
while tf <= tf_max
    % Solve the problem and propagate the solution
    [p0, rho] = solve_fixed_time_transfer_indirect(P, t0, tf);
    
    % If no solution was found continue to the next iteration
    if isnan(p0)
        tf = tf + dtf
        continue
    end

    [~, X] = propagator_MEE_indirect(P, t0, tf, p0(1:6), p0(7), rho);

    % Calculate delta-V
    dv = calc_delta_v(P, X)

    % Record transfer time, delta-V, and initial adjoints
    tf_list = [tf_list, tf];
    dv_list = [dv_list, dv];
    p0_list = [p0_list, p0];

    tf = tf + dtf
end

% Save results to a file
filename = sprintf('optimal_transfers_simple_%s', datestr(now,'mm-dd-yyyy HH-MM'));
save(filename, 'tf_list', 'dv_list', 'p0_list')