function [z, T] = solve_convection_diffusion(N, v, a, b, Q0, alpha0, Tout, T0)
h = 1/N;
z = (0:N)'*h;
alpha = sqrt(v^2/4 + alpha0^2) - v/2;

left = -1/h^2 - v/(2*h);
middle = 2/h^2;
right = -1/h^2 + v/(2*h);

A = sparse(N, N);
rhs = zeros(N, 1);

for j = 1:N - 1
    zj = j*h;

    if j == 1
        rhs(j) = rhs(j) - left*T0;
    else
        A(j, j - 1) = left;
    end

    A(j, j) = middle;
    A(j, j + 1) = right;

    if zj >= a && zj <= b
        rhs(j) = rhs(j) + Q0*sin((zj - a)*pi/(b - a));
    end
end

A(N, N - 2) = 1;
A(N, N - 1) = -4;
A(N, N) = 3 + 2*h*alpha;
rhs(N) = 2*h*alpha*Tout;

U = A\rhs;
T = [T0; U];
end
