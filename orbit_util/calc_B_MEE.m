function B = calc_B_MEE(x, mu)
% Calculates the matrix that maps instantaneous acceleration to time rate
% of change of modified equinoctial elements (MEEs). This is the B matrix
% in the equation
%   x_mee_dot = A + B*u
% 
% Units expected are km, radians, and km^3/s^2.

p = x(1);
f = x(2);
g = x(3);
h = x(4);
k = x(5);
L = x(6);

% Precalculate and cache some things for speed
sL = sin(L);
cL = cos(L);
hsLkcL = h*sin(L) - k*cos(L);

% Intermediate quantitites
w = 1 + f*cL + g*sL;
winv = 1/w;
s2 = 1 + h^2 + k^2;

% Calculate B matrix
B = sqrt(p/mu)*[
    0, 2*p/w, 0;
    sL, winv*((w+1)*cL + f), -(g/w)*hsLkcL;
    -cL, winv*((w+1)*sL + g), (f/w)*hsLkcL;
    0, 0, s2*cL/(2*w);
    0, 0, s2*sL/(2*w);
    0, 0, winv*hsLkcL
];