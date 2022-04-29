function [t, x] = rk4_propagator(f, x0, t0, tf, dt)
% Simple fixed timestep RK4 propagator

t = t0:dt:tf;
N = length(x0);
x = zeros(N, length(t));

for i = 1:length(t)
    k1 = f(t(i), x0);
    k2 = f(t(i) + dt/2, x0 + k1*dt/2);
    k3 = f(t(i) + dt/2, x0 + k2*dt/2);
    k4 = f(t(i) + dt, x0 + k3*dt);
    
    x0 = x0 + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
    x(:,i) = x0;
end