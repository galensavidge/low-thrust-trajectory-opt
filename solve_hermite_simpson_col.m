function [xk, uk, t, x, u] = solve_hermite_simpson_col(f_func, ...
    K_func, L_func, u_max, x0_func, xf_func, tk, xk0, uk0, dt)
% Solves the separated form of the Hermite-Simpson collocation method.
% 
% Arguments:
%   f_func: A function that returns xdot as a function f(x, u), where xdot
%       and x are size (xdim, 1), and u is size (udim, 1).
%   K_func: Terminal cost component. Must be a function of (x0, xf).
%   L_func: Integral cost component. Must be a function of (x, u).
%   u_max: Maximum allowable control norm.
%   x0_func: Equality constraint function for initial state. Should return
%       zero when the constraint is satisfied.
%   xf_func: Equality constraint function for final state.
%   tk: List of sampling times for collocation points.
%   xk0: Column vector of (xdim, 2*N - 1) initial state guesses.
%   uk0: Column vector of (udim, 2*N - 1) initial control guesses.
%   dt: Timestep used for interpolated points

xdim = size(xk0, 1);
udim = size(uk0, 1);

% One-dimensional vector of decisions variables
X0 = [reshape(xk0, [], 1); reshape(uk0, [], 1)];

cost_wrapper = @(X)(evaluate_cost(X, tk, K_func, L_func, xdim, udim));
constraint_wrapper = @(X)(constraint_functions(X, f_func, u_max, ...
    x0_func, xf_func, tk, xdim, udim));

% Optimization problem
problem = struct();
problem.objective = cost_wrapper;
problem.x0 = X0;
problem.lb = [];
problem.ub = [];
problem.nonlcon = constraint_wrapper;
problem.solver = 'fmincon';
options = optimoptions('fmincon', 'Display', 'iter-detailed', ...
    'MaxIterations', 200, 'MaxFunctionEvaluations', Inf, ...
    'UseParallel', true);
problem.options = options;
X_opt = fmincon(problem);

N = length(tk);
M = xdim*(2*N-1);  % Collocation points and "1/2" points
xk = X_opt(1:M);
uk = X_opt(M+1:end);

xk = reshape(xk, xdim, []);
uk = reshape(uk, udim, []);

[t, x, u] = interpolate_collocation_points(tk, xk, uk, f_func, dt);

end

function cost = evaluate_cost(X, tk, K_func, L_func, xdim, udim)
    N = length(tk);
    dt = tk(2:end) - tk(1:end-1);
    
    x0 = vec_from_flat_list(X, 1, xdim);
    xf = vec_from_flat_list(X, 2*N-1, xdim);
    
    % Extract state and control column vectors
    M = xdim*(2*N-1);  % Collocation points + "1/2" points
    xk = X(1:M);
    uk = X(M+1:end);
    
    integral_cost = 0;
    for i=1:N-1
        % Get delta t at i
        dti = dt(i);
        ii = 2*i-1;

        % Get state and control vectors at i, i+1/2, and i+1
        xi = vec_from_flat_list(xk, ii, xdim);
        ui = vec_from_flat_list(uk, ii, udim);
        xi12 = vec_from_flat_list(xk, ii+1, xdim);
        ui12 = vec_from_flat_list(uk, ii+1, udim);
        xi1 = vec_from_flat_list(xk, ii+2, xdim);
        ui1 = vec_from_flat_list(uk, ii+2, udim);
        
        % Evaluate cost function at i, i+1/2, and i+1
        Li = L_func(xi, ui);
        Li12 = L_func(xi12, ui12);
        Li1 = L_func(xi1, ui1);
        
        % Approximate integral
        integral_cost = integral_cost + (dti/6)*(Li + 4*Li12 + Li1);
    end

    cost = K_func(x0, xf) + integral_cost;
end

function [c, ceq] = constraint_functions(X, f_func, u_max, x0_func, ...
    xf_func, tk, xdim, udim)
    % Set up column vectors of inequality and equality constraints
    c = control_constraints(X, u_max, tk, xdim, udim);
    ceq = [dynamics_constraints(X, tk, f_func, xdim, udim);
        terminal_constraints(X, tk, x0_func, xf_func, xdim)];
    if any(imag(c) ~= 0) || any(imag(ceq) ~= 0)
        pause
    end
end

function t_const = terminal_constraints(X, tk, x0_func, xf_func, xdim)
    % Equality constraints on initial and terminal state
    N = length(tk);

    x0 = vec_from_flat_list(X, 1, xdim);
    xf = vec_from_flat_list(X, 2*N-1, xdim);
    
    x0_const = x0_func(x0);
    xf_const = xf_func(xf);

    t_const = [x0_const; xf_const];
end

function u_const = control_constraints(X, unorm_max, tk, xdim, udim)
    % Inequality constraint to enforce |u_k| <= u_max

    N = length(tk);

    % Extract control column vector
    M = xdim*(2*N-1);  % Collocation points + "1/2" points
    uk = X(M+1:end);
    
    u_const = zeros(N, 1);
    for i=1:2*N-1
        ui = vec_from_flat_list(uk, i, udim);
        u_const(i) = norm(ui) - unorm_max;
    end
end

function x_const = dynamics_constraints(X, tk, f, xdim, udim)
    % Path dynamics equality constraints to ensure that the solution 
    % satisfies the state dynamics equations

    N = length(tk);
    dt = tk(2:end) - tk(1:end-1);

    % Extract state and control column vectors
    M = xdim*(2*N-1);  % Collocation points + "1/2" points
    xk = X(1:M);
    uk = X(M+1:end);

    % Loop over all constraints
    x_const = zeros(xdim*(2*N-2), 1);
    for i=1:N-1
        % Get delta t at i
        dti = dt(i);

        % Index for state and control since they include "1/2" points
        ii = 2*i-1;

        % Get state and control vectors at i, i+1/2, and i+1
        xi = vec_from_flat_list(xk, ii, xdim);
        ui = vec_from_flat_list(uk, ii, udim);
        xi12 = vec_from_flat_list(xk, ii+1, xdim);
        ui12 = vec_from_flat_list(uk, ii+1, udim);
        xi1 = vec_from_flat_list(xk, ii+2, xdim);
        ui1 = vec_from_flat_list(uk, ii+2, udim);
        
        % Evaluate dynamics function at i, i+1/2, and i+1
        fi = f(xi, ui);
        fi12 = f(xi12, ui12);
        fi1 = f(xi1, ui1);
        
        % Constraints for x_i and x_i+1/2
        x_const(1+(ii-1)*xdim:ii*xdim) = ...
            xi1 - xi - (dti/6)*(fi + 4*fi12 + fi1);
        x_const(1+ii*xdim:(ii+1)*xdim) = ...
            xi12 - 0.5*(xi + xi1) - 0.125*dti*(fi - fi1);

        if any(imag(x_const) ~= 0)
            pause
        end
    end
end

function [t_points, x_points, u_points] = ...
    interpolate_collocation_points(tk, xk, uk, f, dt)

    N = length(tk);
    hk = tk(2:end) - tk(1:end-1);

    t0 = tk(1);
    tf = tk(end);
    t_points = t0:dt:tf;

    xdim = size(xk, 1);
    udim = size(uk, 1);

    x_points = zeros(xdim, length(t_points));
    u_points = zeros(udim, length(t_points));
    
    % Loop over collocation points
    for i = 1:N-1
        % Get delta t at i
        hi = hk(i);

        % Index for state and control since they include "1/2" points
        ii = 2*i-1;

        % Get state and control vectors at i, i+1/2, and i+1
        xi = xk(:,ii);
        ui = uk(:,ii);
        xi12 = xk(:,ii+1);
        ui12 = uk(:,ii+1);
        xi1 = xk(:,ii+2);
        ui1 = uk(:,ii+2);
        
        % Calculate dynamics at i, i+1/2, and i+1
        fi = f(xi, ui);
        fi12 = f(xi12, ui12);
        fi1 = f(xi1, ui1);
        
        % Time range for this segment
        t_seg = tk(i):dt:(tk(i+1)-dt);
        
        % Calculate state and control over this segment
        x_seg = zeros(xdim, length(t_seg));
        u_seg = zeros(udim, length(t_seg));

        for j = 1:length(t_seg)
            t = t_seg(j);

            tau = t-tk(i);
            
            % Cubic state interpolating function
            x_seg(:,j) = xi + hi*(fi*tau/hi ...
                + 0.5*(-3*fi + 4*fi12 - fi1)*(tau/hi)^2 ...
                + (1/3)*(2*fi - 4*fi12 + 2*fi1)*(tau/hi)^3);
            
            % Quadratic control interpolating function
            u_seg(:,j) = (2/hi^2)*(tau - hi/2)*(tau - hi)*ui ...
                - (4/hi^2)*tau*(tau - hi)*ui12 ...
                + (2/hi^2)*tau*(tau - hi/2)*ui1;
        end

        point_range = 1+(i-1)*length(t_seg):i*length(t_seg);
        x_points(:,point_range) = x_seg;
        u_points(:,point_range) = u_seg;
    end

end

function ri = vec_from_flat_list(rk, idx, dim)
    % Gets vector number "idx" of length "dim" from "rk", a one-dimensional
    % list of flattened vectors
    ri = rk(1+(idx-1)*dim:idx*dim);
end