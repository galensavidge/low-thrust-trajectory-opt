function T = calc_thrust_MEE_indirect(x, m, px, pm, mu, g0, Isp, Tmax, rho)

% Switching function
B = calc_B_MEE(x, mu);
S = (1/m)*norm(B'*px) + pm/(g0*Isp);

% Control smoothing
T = 0.5*Tmax*(1 + tanh(S/rho));