function dv = calc_delta_v(P, X)

% Compute non-dimensionalized version
mf = X(7,end);
dv = P.Isp*P.g0*log(1/mf);

% Convert to km/s
dv = dv*P.LU/P.TU;