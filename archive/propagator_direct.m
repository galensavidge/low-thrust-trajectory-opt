function [t, r, v] = propagator_direct(r0, v0, t0, tf, u, mu)

x0 = [r0; v0];

dynamics_func = @(t, x)(dyanmics_twobody(x, u, mu));

options = odeset('RelTol',1e-9);
[t, x] = ode45(dynamics_func, [t0 tf], x0, options);

r = x(:,1:3)';
v = x(:,4:6)';

end