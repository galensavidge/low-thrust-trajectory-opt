function make_control_plot(P, t, X, rho)

x = X(1:6,:);
m = X(7,:);
px = X(8:13,:);
pm = X(14,:);

[T, u] = control_from_MEE_adjoints(P, x, m, px, pm, rho);

% Convert units
T = T*P.MU*P.LU/(P.TU^2)*1e3;
day = 86400/P.TU;
t = t/day;

figure()
hold on
grid on
plot(t, T)
xlabel('Time [days]', 'Interpreter', 'latex')
ylabel('Thrust [N]', 'Interpreter', 'latex')

figure()
for i = 1:3
    subplot(3, 1, i)
    hold on
    grid on
    plot(t, u(i,:))
end