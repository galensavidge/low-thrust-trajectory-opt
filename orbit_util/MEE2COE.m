function coe = MEE2COE(mee)
% Transforms modified equinoctial elements into classical oribital
% elements. Input angles are in radians. Output angles are in degrees.

p = mee(1);
f = mee(2);
g = mee(3);
h = mee(4);
k = mee(5);
L = mee(6);

% Eccentricity
e = sqrt(f^2 + g^2);

% SMA
a = p/(1-e^2);

% RAAN
RAAN = atan2d(k, h);

% Inclination
i = 2*atand(sqrt(h^2 + k^2));

% AOP
AOP = atan2d(g, f) - RAAN;

% True anomaly
true_anom = rad2deg(L) - RAAN - AOP;

coe = [a; e; i; RAAN; AOP; true_anom];