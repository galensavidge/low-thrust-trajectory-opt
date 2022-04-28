clear
clc
close all
format short g

addpath('..')

% Defined as [velocity; position]
x0 = [1; 1];
xf = [0; 0];

u_max = 1;

% Collocation time points
t0 = 0;
tf = 10;
N = 20;
dt = (tf-t0)/(N-1);
tk = t0:dt:tf;

% Linearly interpolate state for initial guess
xk0 = [x0];
for i=1:N-1
    xki = x0 + (xf - x0)*(dt*i)/(tf - t0);
    xki12 = x0 + (xf - x0)*(dt*i - dt/2)/(tf - t0);
    xk0 = [xk0, xki12, xki];
end

xk0 = zeros(2, 2*N-1);

% Some constant initial guess for control
uk0 = -0.1*ones(1, 2*N-1);

% Wrapper functions for the solver to use
dynamics_func = @(x, u)([u; x(1)]);  % Simple 2nd order system
x0_func = @(x)(x - x0);
xf_func = @(x)(x - xf);
K_func = @(x0, xf)(0);
L_func = @(x, u)(u^2);

[xk, uk, t, x, u] = solve_hermite_simpson_col(dynamics_func, K_func, L_func, ...
    u_max, x0_func, xf_func, tk, xk0, uk0, 0.001);

tk_ext = sort([tk, (tk(1:end-1) + tk(2:end))/2]);

figure()
hold on
scatter(tk_ext, xk(1,:))
scatter(tk_ext, xk(2,:))
plot(t, x(1,:))
plot(t, x(2,:))
legend('$\dot{x}$', '$x$', 'Interpreter', 'latex')

figure()
hold on
scatter(tk_ext, uk)
plot(t, u)