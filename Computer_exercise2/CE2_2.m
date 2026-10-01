clc, clear all, close all 
%Part 2

%a 
T = 40;
Text = 25;
Ly = 5;
Lx = 12;
N = 120;

f = @(x, y) 100 * exp(- 1 / 2 .* (x - 4).^2 - 4 .* (y - 1).^2);

deltaT = T / N;
t = (0:N) * deltaT;

[A, b, x, y, h, M] = Build_system(N, Lx, Ly, Text, f);
u = zeros((N + 1) * M , N + 1);
u_0 = f(x,y(2:end)')';
u(:,1) = Text * ones(length(b),1) + u_0(:);

for i = 1:length(t) - 1 
    u(:, i + 1) = cranknicolson(u(:, i), A, b, deltaT, h, N, M);
end 

%---------------------------------------
%part b 

u_0 = [ones(N + 1, 1) * Text, reshape(u(:, 1), N + 1, [])]; 
u_1 = [ones(N + 1, 1) * Text, reshape(u(:, round(1 / T * N) + 1), N + 1, [])];
u_4 = [ones(N + 1, 1) * Text, reshape(u(:, round(4 / T * N) + 1), N + 1, [])];
u_12 = [ones(N + 1, 1) * Text, reshape(u(:, round(12 / T * N) + 1), N + 1, [])];
u_22 = [ones(N + 1, 1) * Text, reshape(u(:, round(22 / T * N) + 1), N + 1, [])];
u_40 = [ones(N + 1, 1) * Text, reshape(u(:, end), N + 1, [])];

figure(1)
mesh(x, y, u_0')
xlabel('x')
ylabel('y')
title('u(x,y,0)')

figure(2)
mesh(x, y, u_1')
xlabel('x')
ylabel('y')
title('u(x,y,1)')

figure(3);
mesh(x, y, u_4')
xlabel('x')
ylabel('y')
title('u(x,y,4)')

figure(4);
mesh(x, y, u_12')
xlabel('x')
ylabel('y')
title('u(x,y,12)')

figure(5);
mesh(x, y, u_22')
xlabel('x')
ylabel('y')
title('u(x,y,22)')

figure(6);
mesh(x, y, u_40')
xlabel('x')
ylabel('y')
title('u(x,y,40)')
%--------------------------------------
%part c

xidx = 6 / Lx * N + 1;
yidx = 2 / Ly * M + 1;
u_t = u((2 / Ly * M - 1) * (N + 1) + 6 / Lx * N + 1, :);

figure(7)
hold on
plot(t, u_t, 'LineWidth', 1.2)
plot(t, 47.224999298104 * ones(length(t)), '--', 'LineWidth', 1.2)
fprintf('u(6,2,40) = %.12f\n', u_t(end));
legend('u(6,2,t)','T(6,2)', 'Location', 'southeast')

%---------------------------------------
%functions 

function unext = cranknicolson(u, A, b, deltaT, h, N, M)
C = (speye(size(A)) - deltaT / 2 * A);
E = (speye(size(A)) + deltaT / 2 * A) * u + b * deltaT;
unext = C \ E;

end 

function [A, b, x, y, h, M] = Build_system(N, Lx, Ly, Text, f)
h = Lx/N;
M = Ly/h;

x = (0:N)*h;
y = (0:M)*h;

number_unknowns = (N + 1) * M;
A = spalloc(number_unknowns, number_unknowns, 5*number_unknowns);
b = zeros(number_unknowns, 1);

for j = 1:M
    for i = 0:N
        r = i + 1 + (j - 1) * (N + 1);
        A(r, r)= - 4 / h^2;
        
        if i == 0 %Neoman condition
            A(r, r + 1) = 2 / h^2;
        elseif i == N %Neoman condition
            A(r, r - 1) = 2 / h^2; 
        else
            A(r, r - 1) = 1 / h^2; 
            A(r, r + 1) = 1 / h^2;
        end

        if j == 1 % insulated lower boundrary
            b(r) = Text / h^2;
        else
            A(r, r -(N + 1)) = 1 / h^2;
        end

        if j == M % insulated upper boundrary
            A(r, r - (N + 1)) = A(r, r - (N + 1)) + 1 / h^2;
        else
            A(r, r + (N + 1)) = 1 / h^2;
        end
        b(r) = b(r) + f(i * h, j * h); % f(x_i, y_i)
    end
end

end 