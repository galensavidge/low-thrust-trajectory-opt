function plot_MEE_as_cartesian(t, x)

r =  position_from_MEE(x);

hold on
grid on
axis equal
plot3(r(1,:), r(2,:), r(3,:))