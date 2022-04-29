function uhat = calc_control_dir_MEE_indirect(x, px, mu)

% Find control direction from primer vector
B = calc_B_MEE(x, mu);
v = B'*px;
uhat = -v/norm(v);