function [T, u] = control_from_MEE_adjoints(P, x, m, px, pm, rho)
% Calculates the LVLH control vector from state and adjoint history

u = zeros(3, length(x));
T = zeros(1, length(x));
for i = 1:length(u)
    T(i) = calc_thrust_MEE_indirect(x(:,i), m(i), px(:,i), pm(i), P.mu, ...
        P.g0, P.Isp, P.Tmax, rho);
    uhat = calc_control_dir_MEE_indirect(x(:,i), px(:,i), P.mu);
    u(:,i) = T(i)*uhat;
end