function make_pareto_front_plot(P, tf_list, dv_list)

% Convert units
day = 86400/P.TU;
tf_list = tf_list/day;

figure()
hold on
grid on
scatter(tf_list, dv_list)
xlabel('Transfer Time [days]', 'Interpreter', 'latex')
ylabel('$\Delta V [km/s]$', 'Interpreter', 'latex')