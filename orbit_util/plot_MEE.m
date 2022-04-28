function plot_MEE(t, x)

title('Modified Equinoctial Orbital Elements')


subplot(2, 3, 1)
hold on
plot(t, x(1,:))
grid on
title('p')

subplot(2, 3, 2)
hold on
plot(t, x(2,:))
grid on
title('f')

subplot(2, 3, 3)
hold on
plot(t, x(3,:))
grid on
title('g')

subplot(2, 3, 4)
hold on
plot(t, x(4,:))
grid on
title('h')

subplot(2, 3, 5)
hold on
plot(t, x(5,:))
grid on
title('k')

subplot(2, 3, 6)
hold on
plot(t, x(6,:))
grid on
title('L')