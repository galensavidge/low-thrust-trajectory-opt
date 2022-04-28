function P = problem_setup()

% Initial GTO
x0 = [11363.195; 0.73136656; 0; 0.25396765; 0; 1.5707963];

% Final GEO
xf = [42164; 0; 0; 0; 0; 0];  %  Note Lf is unconstrained

% Initial mass
m0 = 370;  % kg

% Earth gravitation parameter
mu = 3.986004418e5;  % km^3/s^2

% EP thruster parameters
g0 = 9.80665e-3;  % km/s^2
Isp = 1640;  % sec
Tmax = 68e-6;  % kg-km/s^2

% Normalize things to the range [0, 1]
LU = xf(1);  % Measure length in units of GEO radius
TU = LU/sqrt(mu/LU);  % Set up time units such that mu = 1
MU = m0;  % Measure mass in units of spacecraft initial mass

mu = mu*TU^2/LU^3;  % Will always be 1, but good sanity check

% Normalize distances
x0(1) = x0(1)/LU;
xf(1) = xf(1)/LU;

% Normalize mass
m0 = m0/MU;

% Normalize thruster parameters
g0 = g0*TU^2/LU;
Isp = Isp/TU;
Tmax = Tmax*TU^2/(MU*LU);

% Store things in a struct and return
P = struct();
P.x0 = x0;
P.xf = xf;
P.m0 = m0;
P.mu = mu;
P.g0 = g0;
P.Isp = Isp;
P.Tmax = Tmax;
P.LU = LU;
P.TU = TU;
P.MU = MU;