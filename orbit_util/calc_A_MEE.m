function A = calc_A_MEE(x, mu)
% Calculates the A matrix in
%   x_mee_dot = A + B*u
% 
% Units expected are km, radians, and km^3/s^2.

p = x(1);
f = x(2);
g = x(3);
L = x(6);

w = 1 + f*cos(L) + g*sin(L);
A = [0; 0; 0; 0; 0; sqrt(mu*p)*(w/p)^2];