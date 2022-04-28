clear
clc
close all
format long g

addpath('..')

% Problem setup
% Defined as [velocity; position]
x0 = [1; 1];
xf = [0; 0];

u_max = 1;

t0 = 0;
tf = 10;

p0_guess = [1; -2];  % [Pv; Pu]

f = @(p0)(single_shot(x0, p0, t0, tf, xf));

p0_guess = random_sample_function(f, [-1; 0], [1; 1], 1000);
p0_soln = newton_iteration(f, p0_guess, 0.1, 1e-3, 1000, 1e-9)
[t, x, p, u] = propagate_dynamics(x0, p0_soln, t0, tf);

figure()
hold on
grid on
plot(t, x(1,:))
plot(t, x(2,:))
title('States')
legend('$x$', '$\dot{x}$', Interpreter='latex')

figure()
hold on
grid on
plot(t, p(1,:))
plot(t, p(2,:))
title('Adjoints')
legend('$p_v$', '$p_u$', Interpreter='latex')

figure()
hold on
grid on
plot(t, u)
title('Control')

% Piecewise function for optimal control in terms of adjoints
function u = u_opt(p)
    pu = p(2);
    if abs(pu) < 1
        u = 0;
    else
        u = -sign(pu);
    end
end

function u = u_opt_energy(p)
    u = -0.5*p(2);
end

function xdot = state_dynamics(x, p)
    u = u_opt_energy(p);
    xdot = [x(2); u];
end

function pdot = adjoint_dynamics(p)
    pdot = [0; p(1)];
end

function [t, x_hist, p_hist, u_hist] = propagate_dynamics(x0, p0, t0, tf)
    % Vector including state and adjoints
    X0 = [x0; p0];

    Xdot_func = @(t, X)([state_dynamics(X(1:2), X(3:4)); ...
        adjoint_dynamics(X(3:4))]);

    [t, X_hist] = ode45(Xdot_func, [t0 tf], X0);

    t = t';
    X_hist = X_hist';

    x_hist = X_hist(1:2,:);
    p_hist = X_hist(3:4,:);
    u_hist = zeros(1, length(t));
    for i=1:length(t)
        u_hist(i) = u_opt_energy(p_hist(:,i));
    end
end

function xf_err = single_shot(x0, p0, t0, tf, xf)
    [~, x_hist, ~, ~] = propagate_dynamics(x0, p0, t0, tf);
    xf_err = x_hist(:,end) - xf;
end