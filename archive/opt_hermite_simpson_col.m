clear
clc
close all
format long e

% mu = 3.986004418e5;

% R_earth = 6378.1;  % km

% Initial GTO
% x0 = [11363.195; 0.73136656; 0; 0.25396765; 0; 1.5707963];
% x0 = [11363.195; 0; 0; 0; 0; 1.5707963];
% Final GEO
% xf = [42164; 0; 0; 0; 0; 3*pi/2];  %  Note L is unconstrained

[x0, m0, xf, mu, g0, Isp, Tmax, TU] = problem_setup();

% Constrain transfer time
t0 = 0;
tf = 90*86400/TU;  % sec

% Set up evenly spaced collocation points
num_points = 30;
dt = (tf-t0)/(num_points-1);
tk = t0:dt:tf;

% Linearly interpolate state for initial guess
xk0 = [x0];
for i=1:num_points-1
    xki = x0 + (xf - x0)*(dt*i)/(tf - t0);
    xki12 = x0 + (xf - x0)*(dt*i - dt/2)/(tf - t0);
    xk0 = [xk0, xki12, xki];
end

% Spacecraft properties
M_sc = 370;  % kg
F_BOL = 68e-3;  % Beginning of life max thrust, N-m
u_max = F_BOL / M_sc * 1e-3;  % km/s^2

uk0 = zeros(3, 2*num_points - 1);

dynamics_func = @(x, u)(dynamics_MEE(x, u, mu));

cost_func = @(xk, uk, tk)(u_squared_cost(uk, u_max));

x0_func = @(x)(initial_constraint(x, x0));
xf_func = @(x)(terminal_constraint(x, xf));

L_func = @(x, u)(sum(u.^2));
K_func = @(x, u)(0);
dt = 1;  % sec

[xk, uk, t, x, u] = solve_hermite_simpson_col(dynamics_func, K_func, L_func, u_max, ...
    x0_func, xf_func, tk, xk0, uk0, dt);

tk_ext = sort([tk, (tk(1:end-1) + tk(2:end))/2]);
plot_MEE(tk_ext, xk)

plot_MEE(t, x)

r = zeros(3, length(t));
for i = 1:length(t)
    r(:,i) = position_from_MEE(x(:,i));
end

plot_position(t, r)

function cost = initial_constraint(x, x0)
    % Constrain initial state to x0
    cost = [(x(1) - x0(1))/x0(1); x(2:6) - x0(2:6)];
end

function cost = terminal_constraint(x, xf)
    % Constrain final state to xf, except L
    cost = [(x(1) - xf(1))/xf(1); x(2:5) - xf(2:5)];
end