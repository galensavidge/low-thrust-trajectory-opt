function [t, x] = propagator_MEE_direct(x0, t0, tf, u, mu)

dynamics_func = @(t, x)dynamics_MEE(x, u, mu);

options = odeset('RelTol',1e-8,'AbsTol',1e-9);
[t, x] = ode45(dynamics_func, [t0 tf], x0, options);
t = t';
x = x';