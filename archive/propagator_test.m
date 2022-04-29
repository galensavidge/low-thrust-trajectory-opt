clear
clc
close all

addpath('orbit_util')

mu = 3.986004418e5;

R_earth = 6378.1;  % km

% Circular LEO
x0 = [R_earth + 500; 0; 0; 0; 0; 0];

% Constrain transfer time
tf = 15000;  % [sec]

num_segments = 3;
segment_time = tf/num_segments;
segment_times = 0:segment_time:tf;
u = [0; 0; 0; 0; 0; 0.003; 0; 0; 0];

[t, x] = propagator_MEE_thrust_segments(x0, segment_times, u, mu);

% Calculate cartesian position from MEEs
r = zeros(3, length(t));
for j = 1:length(t)
    r(:,j) = position_from_MEE(x(:,j));
end

plot_MEE(t, x)
plot_position(t, r)