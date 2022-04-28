clear
clc
close all

mu = 3.986004418e5;

R_earth = 6378.1;  % km

% Circular orbit
x0 = [R_earth + 500; 0.5; 0; 0; 0; 0];

% Constrain transfer time
tf = 86400;  % [sec]

rf = R_earth + 2500;

num_segments = 100;
segment_time = tf/num_segments;
segment_times = 0:segment_time:tf;
% Constant thrust arcs in the LVLH frame
u = 1e-4*ones(3*num_segments, 1);
u_max = 1e-3;

% Cost and constraint functions
cost_func = @(u)(calc_cost(u, segment_times));
constraint_func = @(u)(nl_constraint(x0, segment_times, u, mu, rf));

% Set up optimizer
problem = struct();
problem.objective = cost_func;
problem.x0 = u;
problem.lb = zeros(3*num_segments, 1);
problem.ub = u_max*ones(3*num_segments, 1);
problem.nonlcon = constraint_func;
problem.solver = 'fmincon';
options = optimoptions('fmincon', 'Display', 'iter-detailed', ...
    'MaxIterations', 100, 'MaxFunctionEvaluations', Inf, ...
    'UseParallel', true);
problem.options = options;
uopt = fmincon(problem);

[t, x] = propagator_MEE_thrust_segments(x0, segment_times, uopt, mu);

% Calculate cartesian position from MEEs
r = zeros(3, length(t));
for j = 1:length(t)
    r(:,j) = position_from_MEE(x(:,j));
end

plot_MEE(t, x)
plot_position(t, r)

function dv = calc_cost(u, segment_times)
    % Simple integration of norm of constant thrust segments
    u_norms = sqrt(u(1:3:end).^2 + u(2:3:end).^2 + u(3:3:end).^2)';
    dv = sum(u_norms.*(segment_times(2:end) - segment_times(1:end-1)));
end

function [c, ceq] = nl_constraint(x0, segment_times, u_list, mu, rf)
    % Unused inequality constraint
    c = [0];

    % Propagate to tf
    [~, x] = propagator_MEE_thrust_segments(x0, segment_times, u_list, mu);
    xf = x(:,end);

    % Extract MEEs
    pf = xf(1);
    ff = xf(2);
    gf = xf(3);
    hf = xf(4);
    kf = xf(5);
    
    % We want to drive p to rf, and the others to zero for a circular
    % orbit. Lf is unconstrained.
    ceq = [(pf - rf)/rf; ff; gf; hf; kf];
end
