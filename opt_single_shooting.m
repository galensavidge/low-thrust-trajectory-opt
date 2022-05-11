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
tf = 180*day;  % Fixed transfer time

% Solve 2PBVP for optimal transfer
[best_p0, rho] = solve_fixed_time_transfer_indirect(P, t0, tf);
% best_p0 = [-0.320261822914485;-0.193738344243347;6.555566905522959e-05;
%     0.235110780657784;-1.600528169647235e-04;-7.297911479765603e-06;
%     -0.787678616604158];
% rho = 1e-5;

% Propagate solution
[t, X] = propagator_MEE_indirect(P, t0, tf, best_p0(1:6), best_p0(7), rho);

x = X(1:6,:);
m = X(7,:);
px = X(8:13,:);
pm = X(14,:);

% Plot control
make_control_plot(P, t, X, rho)

% Plot equinoctial elements over time
make_MEE_plot(P, t, X)

% Plot cartesian transfer trajectory
make_cartesian_plot(P, X, rho)

% Find final errors
xf = x(:, end);
xf_err = xf - P.xf;
xf_err(1) = xf_err(1)*P.LU;

px_err = px(6,end);
pm_err = pm(end) + 1;