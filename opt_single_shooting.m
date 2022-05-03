clear
clc
close all
format long g

addpath('orbit_util')
addpath('orbit_util\autogen\')

% Set up the problem
P = problem_setup();

t0 = 0;
day = 86400;  % One day in seconds
tf = 180*day/P.TU;

% [best_p0, rho] = solve_fixed_time_transfer_indirect(P, t0, tf);

best_p0 = [-0.320261822914485;-0.193738344243347;6.555566905522959e-05;
    0.235110780657784;-1.600528169647235e-04;-7.297911479765603e-06;
    -0.787678616604158];
rho = 1e-5;


% Generate Pareto front
% tf_min = 120*day/P.TU;
% tf_max = 240*day/P.TU;
% dtf = 15*day/P.TU;
% [tf_list, dv_list, p0_list] = build_pareto_front(P, t0, tf, rho, ...
%     best_p0, tf_min, tf_max, dtf);
% 
% figure()
% scatter(tf_list, dv_list)

filename = sprintf('optimal_transfers_%s', datestr(now,'mm-dd-yyyy HH-MM'));
save(filename, 'tf_list', 'dv_list', 'p0_list')

% Plot out solution
[t, X] = propagator_MEE_indirect(P, t0, tf, best_p0(1:6), best_p0(7), rho);

x = X(1:6,:);
m = X(7,:);
px = X(8:13,:);
pm = X(14,:);

figure()
plot_MEE(t, x);

[T, u] = control_from_MEE_adjoints(P, x, m, px, pm, rho);

figure()
hold on
grid on
plot(t, T)

figure()
for i = 1:3
    subplot(3, 1, i)
    hold on
    grid on
    plot(t, u(i,:))
end

r = position_from_MEE(x);
plot_position(t, r);
