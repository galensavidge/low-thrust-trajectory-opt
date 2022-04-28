function plot_position(t, r)
% Plots keplerian position. Time history is of length N while r is 3xN.

figure()
hold on
grid on
axis equal
plot3(r(1,:), r(2,:), r(3,:))

figure()
subplot(1, 3, 1)
hold on
grid on
plot(t, r(1,:))
subplot(1, 3, 2)
hold on
grid on
plot(t, r(2,:))
subplot(1, 3, 3)
hold on
grid on
plot(t, r(3,:))