function [t, y] = impeuler(Tend, y0, h)
%IMPEULER Solve Robertson problem with the Implicit Euler method.

r1 = 0.05;
r2 = 1.2e4;
r3 = 4e7;

f = @(u) [-r1*u(1) + r2*u(2)*u(3);
           r1*u(1) - r2*u(2)*u(3) - r3*u(2)^2;
           r3*u(2)^2];

Jf = @(u) [-r1,  r2*u(3),                 r2*u(2);
            r1, -r2*u(3) - 2*r3*u(2),   -r2*u(2);
             0,  2*r3*u(2),                       0];

I = eye(3);
tol = 1e-10;

t = 0:h:Tend;
u = y0(:);
y = u';

for t1 = t(2:end)
    d = 1;
    v = u;

    while norm(d) > tol
        F = v - h*f(v) - u;
        J = I - h*Jf(v);
        d = -J\F;
        v = v + d;
    end

    u = v;
    y = [y; u'];
end

t = t';
end
