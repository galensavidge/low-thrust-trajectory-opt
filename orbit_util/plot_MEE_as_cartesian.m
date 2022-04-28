function plot_MEE_as_cartesian(t, x)

r = zeros(3, length(x));
for i = 1:length(x)
    r(:,i) = position_from_MEE(x(:,i));
end

hold on
grid on
axis equal
plot3(r(1,:), r(2,:), r(3,:))