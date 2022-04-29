function xdot = dyanmics_twobody(x, u, mu)
% Evaluates two-body cartesian equations of motion.

r = x(1:3);
v = x(4:6);

rdot = v;
vdot = -mu/(norm(r)^3)*r + u;

xdot = [rdot; vdot];