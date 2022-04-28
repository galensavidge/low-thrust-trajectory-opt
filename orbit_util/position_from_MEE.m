function r_N = position_from_MEE(x)

% Translate MEE to classical orbital elements (COE)
coe = MEE2COE(x);

% Extract orbital elements from state vectors
p = x(1);
e = coe(2);
i = coe(3);
RAAN = coe(4);
AOP = coe(5);
true_anom = coe(6);

% Construct position in the local orbital (LVLH) frame
r = p/(1 + e*cosd(true_anom));
r_LVLH = [r; 0; 0];

% Express in inertial frame
C = C_LVLH_to_inertial(RAAN, i, AOP, true_anom);
r_N = C*r_LVLH;
