function [x, y, T] = solve_heat_2d(N, Lx, Ly, Text, f)
h = Lx/N;
M = Ly/h;

x = (0:N)*h;
y = (0:M)*h;

number_unknowns = (N + 1)*M;
A = spalloc(number_unknowns, number_unknowns, 5*number_unknowns);
rhs = zeros(number_unknowns, 1);

%% Interior points
for j = 1:M - 1
    for i = 1:N - 1
        k = i + 1 + (j - 1)*(N + 1);

        A(k, k) = 4;
        A(k, k - 1) = -1;
        A(k, k + 1) = -1;
        A(k, k + N + 1) = -1;

        rhs(k) = h^2*f(i*h, j*h);

        if j == 1
            rhs(k) = rhs(k) + Text;
        else
            A(k, k - N - 1) = -1;
        end
    end
end

%% Left and right boundaries
for j = 1:M
    k_left = 1 + (j - 1)*(N + 1);
    A(k_left, k_left) = -3;
    A(k_left, k_left + 1) = 4;
    A(k_left, k_left + 2) = -1;

    k_right = N + 1 + (j - 1)*(N + 1);
    A(k_right, k_right) = 3;
    A(k_right, k_right - 1) = -4;
    A(k_right, k_right - 2) = 1;
end

%% Top boundary
for i = 1:N - 1
    k = i + 1 + (M - 1)*(N + 1);
    A(k, k) = 3;
    A(k, k - N - 1) = -4;
    A(k, k - 2*(N + 1)) = 1;
end

U = A\rhs;

T = zeros(M + 1, N + 1);
T(1, :) = Text;

for j = 1:M
    for i = 0:N
        k = i + 1 + (j - 1)*(N + 1);
        T(j + 1, i + 1) = U(k);
    end
end
end
