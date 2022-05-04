function make_cartesian_plot(P, X, rho)

% Unpack the full state and adjoint vector
x = X(1:6,:);
m = X(7,:);
px = X(8:13,:);
pm = X(14,:);

% Calculate position and control over time
r = position_from_MEE(x)*P.LU;
[T, u] = control_from_MEE_adjoints(P, x, m, px, pm, rho);

% Split up position points into thrusting and coast segments
thrust_on = T>0;
r_thrust = r;
r_thrust(:,~thrust_on) = NaN;
r_coast = r;
r_coast(:,thrust_on) = NaN;

% Nice colors
% https://fanwangecon.github.io/M4Econ/graph/tools/fs_color.html
blue = [114 147 203]./255;
red = [211 94 96]./255;
black = [128 133 133]./255;
green = [132 186 91]./255;
brown = [171 104 87]./255;
purple = [144 103 167]./255;

figure()
hold on
grid on
axis equal
plot3(r_coast(1,:), r_coast(2,:), r_coast(3,:), 'Color', black)
plot3(r_thrust(1,:), r_thrust(2,:), r_thrust(3,:), 'Color', blue)

% Plot the initial and final orbits
x0 = x(:,1);
x0(1) = x0(1)*P.LU;
xf = x(:,end);
xf(1) = xf(1)*P.LU;
plot_MEE_orbit(x0, purple)
plot_MEE_orbit(xf, green)

legend('Coast Arcs', 'Thrust Arcs', 'Initial Orbit', 'Final Orbit', Interpreter='latex')
xlabel('$X$ [km]', 'Interpreter', 'latex')
ylabel('$Y$ [km]', 'Interpreter', 'latex')
zlabel('$Z$ [km]', 'Interpreter', 'latex')