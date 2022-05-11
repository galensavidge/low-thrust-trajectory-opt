function plot_MEE_orbit(x, color)

% Get classical orbital elements
coe = MEE2COE(x);
a = coe(1);
e = coe(2);
i = coe(3);
RAAN = coe(4);
AOP = coe(5);

p = a*(1- e^2);

% Calculate cartesian coordinates for one full orbit
true_anom_list = 0:1:360;
r_points = zeros(3, length(true_anom_list));
for ii = 1:length(true_anom_list)
    true_anom = true_anom_list(ii);
    % Find r vector in rotating frame
    r_rth = [p / (1 + e*cosd(true_anom)), 0, 0];
    % Rotate to inertial frame
    C = C_LVLH_to_inertial(RAAN, i, AOP, true_anom);
    r_points(:,ii) = C*r_rth';
end

% Plot orbit
plot3(r_points(1,:), r_points(2,:), r_points(3,:), 'Color', color, ...
    'LineWidth', 2)