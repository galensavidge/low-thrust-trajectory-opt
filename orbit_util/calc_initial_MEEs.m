clear
clc
close all

R_earth = 6378.1363;  % km

% Define elements here
r_a = 42300;  % km
r_p = R_earth + 185  % km
i = 28.5;  % deg
AOP = 0;
RAAN = 0;
true_anom = 90; % deg; puts t0 at top of the SLR

a = 0.5*(r_a + r_p);
e = (r_a - r_p)/(r_a + r_p);

mee = COE2MEE(a, e, i, AOP, RAAN, true_anom)