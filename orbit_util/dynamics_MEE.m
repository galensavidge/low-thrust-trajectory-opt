function xdot = dynamics_MEE(x, u, mu)
% Returns the rate of change of the modified equinoctial elements (MEE) as
% a function of current elements and LVLH acceleration vector u.

A = calc_A_MEE(x, mu);
B = calc_B_MEE(x, mu);
xdot = A + B*u;