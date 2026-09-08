function u = impeuler_nosave(Tend, y0, h)
r1 = 0.05;
r2 = 1.2e4;
r3 = 4e7;

f = @(u) [-r1*u(1) + r2*u(2)*u(3);
           r1*u(1) - r2*u(2)*u(3) - r3*u(2)^2;
           r3*u(2)^2];

Jf = @(u) [-r1,  r2*u(3), r2*u(2);
            r1, -r2*u(3) - 2*r3*u(2),   -r2*u(2);
             0,  2*r3*u(2), 0];

I = eye(3);
tol = 1e-10;
N = round(Tend/h);
u = y0(:);

for n = 1:N
    d = 1;
    v = u;

    while norm(d) > tol
        F = v - h*f(v) - u;
        J = I - h*Jf(v);
        d = -J\F;
        v = v + d;
    end

    u = v;
end
end
