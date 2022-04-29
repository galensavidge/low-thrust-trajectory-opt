function [t, X_hist] = propagator_MEE_indirect(P, t0, tf, px0, pm0, rho)
% Construct state vector
X0 = [P.x0; P.m0; px0; pm0];

dynamics_func = @(t, X)(dynamics_MEE_indirect(X, P.mu, P.g0, P.Isp, P.Tmax, rho));

% Propagate forward to tf
propagator_options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
[t, X_hist] = ode45(dynamics_func, [t0 tf], X0, propagator_options);
X_hist = X_hist';