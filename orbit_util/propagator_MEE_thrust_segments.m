function [t_hist, x_hist] = propagator_MEE_thrust_segments(x0, segment_times, u_list, mu)
% Propagates modified equinoctial elements over a series of constant thrust
% segments.

t_hist = [];
x_hist = [];

x = x0;
for i=1:(length(segment_times)-1)
    % Pull thrust vector and time boundaries for this segment
    u = u_list(3*i-2:3*i);
    t0 = segment_times(i);
    tf = segment_times(i+1);
     
    % Run propagation
    [t_segment, x_segment] = propagator_MEE_direct(x, t0, tf, u, mu);
    x = x_segment(:,end);

    % Save history
    t_hist = [t_hist, t_segment];
    x_hist = [x_hist, x_segment];
end