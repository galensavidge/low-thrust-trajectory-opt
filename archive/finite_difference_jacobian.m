function dfdx = finite_difference_jacobian(f, x, delta)
% Calculates the jacobian of a function using finite differences.
%
% Args:
%   f: Function of x to differentiate
%   x: Point to differentiate at
%   delta: Magnitude of the finite difference to use

xdim = length(x);
y0 = f(x);
ydim = length(y0);
dfdx = zeros(xdim, ydim);
for i=1:xdim
    x1 = x;
    x1(i) = x1(i) + delta;
    y1 = f(x1);
    dfdx(i,:) = (y1 - y0)/delta;
end