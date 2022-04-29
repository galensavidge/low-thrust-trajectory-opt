function mee = COE2MEE(a, e, i, AOP, RAAN, true_anom)
% Converts classical orbital elements to modified equinoctial elements. 
% Input angles are in and degrees. Output angles are in km and radians.

p = a*(1-e^2);
f = e*cosd(AOP + RAAN);
g = e*sind(AOP + RAAN);
h = tand(i/2)*cosd(RAAN);
k = tand(i/2)*sind(RAAN);
L = deg2rad(AOP + RAAN + true_anom);

mee = [p; f; g; h; k; L];