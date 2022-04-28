function C_N_LVLH = C_LVLH_to_inertial(raan, inc, arg_peri, true_anom)
% Finds rotating LVLH frame to inertial frame rotation matrix. Arguments 
% are assumed to be in degrees.

Theta = arg_peri + true_anom;

C1 = [cosd(raan), sind(raan), 0; -sind(raan), cosd(raan), 0; 0, 0, 1];
C2 = [1, 0, 0; 0, cosd(inc), sind(inc); 0, -sind(inc), cosd(inc)];
C3 = [cosd(Theta), sind(Theta), 0; -sind(Theta), cosd(Theta), 0; 0, 0, 1];

C = C3*C2*C1;
C_N_LVLH = C';