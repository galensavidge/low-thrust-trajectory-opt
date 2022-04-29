function T = calc_thrust_MEE_indirect(x, m, px, pm, mu, g0, Isp, Tmax, rho)
% Calculates optimal thrust magnitude for the indirect, modified 
% equinoctial element based primer vector formulation of the low thrust
% transfer problem.

% Primer vector
B = calc_B_MEE(x, mu);
v = B'*px;

% Switching function
S = (1/m)*norm(v) + pm/(g0*Isp);

% Control smoothing
T = 0.5*Tmax*(1 + tanh(S/rho));