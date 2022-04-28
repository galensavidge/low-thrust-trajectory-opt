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
rho = 1;  % Smoothing parameter

% Random guess
% p0_guess = random_sample_function(shot_wrapper, -50*ones(7, 1), 50*ones(7, 1), 1000)

% Particle swarm
% particleswarm_shot_wrapper = @(p0)(sum(abs(single_shot(P, t0, tf, p0(1:6)', p0(7), rho))));
% output_wrapper = @(p0, state)(ps_output_function(P, p0, state, t0, tf, rho));
% 
% options = optimoptions('particleswarm', Display='iter', ...
%     UseParallel=true, OutputFcn=output_wrapper, FunValCheck='on', ...
%     MaxIterations=50, MaxStallIterations=1000);
% best_p0 = particleswarm(particleswarm_shot_wrapper, 7, -1*ones(1, 7), 1*ones(1, 7), options);
% % best_p0 = [-1.003149272005201e+05,-3.611387510559840e+04,27.894598527105934,2.241824990393914e+04,-36.415324566772100,-4.513841538276425,1.154516057453528e+04]';
% best_p0 = best_p0';

% Randomly guess costates until one converges
max_guesses = 100;
iter_tolerance = 1e-6;
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
    rho = rho*0.1
    [best_p0, xf_err] = solve_2pbvp(P, t0, tf, rho, best_p0);
end

% Plot out solution
[t, X] = indirect_propagate(P, t0, tf, best_p0(1:6), best_p0(7), rho);

x = X(1:6,:);
m = X(7,:);
px = X(8:13,:);
pm = X(14,:);

figure()
plot_MEE(t, x);

r = zeros(3, length(x));
for i = 1:length(x)
    r(:,i) = position_from_MEE(x(:,i));
end

plot_position(t, r);

function [t, X_hist] = ...
    indirect_propagate(P, t0, tf, px0, pm0, rho)
    % Construct state vector
    X0 = [P.x0; P.m0; px0; pm0];

    dynamics_func = @(t, X)(dynamics_MEE_indirect(X, P.mu, P.g0, P.Isp, P.Tmax, rho));

    % Propagate forward to tf
    propagator_options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    [t, X_hist] = ode45(dynamics_func, [t0 tf], X0, propagator_options);
    X_hist = X_hist';

    % Fixed time step propagation
%     dt = 0.005;
%     [t, X_hist] = rk4_propagator(dynamics_func, X0, t0, tf, dt);
end

function xf_constraint = single_shot(P, t0, tf, px0, pm0, rho)
    [~, X_hist] = indirect_propagate(P, t0, tf, px0, pm0, rho);
    Xf = X_hist(:,end);
    
    % Constrain final orbit state = xf, pL_f = 0, and pm_f = -1
    xf_constraint = [Xf(1:5) - P.xf(1:5); Xf(13); Xf(14) + 1];
end

function [best_p0, xf_err] = solve_2pbvp(P, t0, tf, rho, p0_guess)
    shot_wrapper = @(p0)(single_shot(P, t0, tf, p0(1:6), p0(7), rho));

    options = optimset('FinDiffRelStep', 1e-9, 'Display', 'iter-detailed', ...
        'TolX', 1e-15, 'UseParallel', true);
    [best_p0, xf_err] = fsolve(shot_wrapper, p0_guess, options);
end

% Output function for particle swarm. Plots the best solution every
% iteration.
function stop = ps_output_function(P, optimValues, state, t0, tf, rho)
    p0 = optimValues.bestx';

    % Propagate the solution
    [xf, t, X] = indirect_propagate(P, t0, tf, p0(1:6), p0(7), rho);
    x = X(1:6,:);
    m = X(7,:);
    px = X(8:13,:);
    pm = X(14,:);

    % Plot the orbit over time
    figure(1)
    clf
    plot_MEE(t, x)
    figure(2)
    clf
    plot_MEE_as_cartesian(t, x)
    stop = false;

    % Plot the control thrust over time
    figure(3)
    T = zeros(1, length(t));
    for i = 1:length(t)
        T(i) = calc_thrust_MEE_indirect(x(:,i), m(i), px(:,i), pm(i), ...
            P.mu, P.g0, P.Isp, P.Tmax, rho)/P.Tmax;
    end
    plot(t, T)

    figure(4)
    hold on
    grid on
    for i=1:length(p0)
        scatter(optimValues.iteration, p0(i))
    end
    legend(['p0_1', 'p0_2', 'p0_3', 'p0_4', 'p0_5', 'p0_6', 'p0_7'])
end