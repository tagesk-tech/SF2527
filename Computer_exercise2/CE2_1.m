clc, clear all, close all 
%Part 1

%a
T = 2; 
d = 0.35; 
a = 1.2;
Nx = 50;
stabelconstant = 0.99; %under 1 for stablity
tau1 = 1.1;

deltax = 1 / Nx;
deltat = stabelconstant * 1 / (2 * d) * deltax^2;

Nt = round(T / deltat);
t = linspace(0, T, Nt);
x = (0:Nx) * deltax;

[A, b, u, g] = Build_system(Nx, Nt, deltax, d, a);
u(:, 1) = g;

for n = 1:Nt - 1
    u(:, n + 1) = expliciteuler(t(n), u(:,n), A, b(:,n), a, d, deltat, deltax);
end  

u_0_t = boundry(Nt, a, t);
u = [u_0_t; u];

disp('Theoretical stability condition')
delta = char(916); 
tau = char(964);
fprintf(['f(',delta, tau, ',', delta, ' x) = ', delta, tau, '/', delta, ' x^2 = %.12f\n'], 1 / (2 * d));

figure(1)
surf(t, x, u, 'EdgeColor', 'none')
ylabel('x')
xlabel('\tau')

t1 = round((tau1 / T) * Nt + 1);
figure(2)
plot(x, u(:, t1))
xlabel('x')
ylabel('u(x,\tau_1)')
ylim([0, 0.5])

figure(3)
hold on 
plot(t, u(1, :))
plot(t, u(end, :))
xlabel('\tau')
ylabel('u')
legend('u(0,\tau)', 'u(1,\tau)')

%--------------------------------------
%b

N = [100, 200, 400]';
tspan = [0, T];
options1 = odeset('RelTol', 1e-5);

ODE23_time = zeros(length(N), 1);
ODE23s_time = zeros(length(N), 1);
ODE23sJ_time = zeros(length(N), 1);

ODE23_iter = zeros(length(N), 1);
ODE23s_iter = zeros(length(N), 1);
ODE23sJ_iter = zeros(length(N), 1);

iterationvec = zeros(3, 3);

for i = 1:length(N)
    Nx = N(i);
    deltax = 1/Nx; 

    [A, b, u, g] = Build_system(Nx, Nt, deltax, d, a);

    options2 = odeset('RelTol', 1e-5, 'Jacobian', A);

    t_start1 = tic;
    [t23, u23] = ode23(@(t,u) odefcn(t, u, A, b(:,1), a, d, deltax), tspan, g, options1);
    ODE23_time(i) = toc(t_start1);
    u23 = [boundry(length(t23), a, t23); u23'];
    ODE23_iter(i) = length(t23) - 1;

    t_start2 = tic;
    [t23s, u23s] = ode23s(@(t,u) odefcn(t, u, A, b(:,1), a, d, deltax), tspan, g, options1);
    ODE23s_time(i) = toc(t_start2);
    u23s = [boundry(length(t23s), a, t23s); u23s'];
    ODE23s_iter(i) = length(t23s) - 1;

    t_start3 = tic;
    [t23sJ, u23sJ] = ode23s(@(t,u) odefcn(t, u, A, b(:,1), a, d, deltax), tspan, g, options2);
    ODE23sJ_time(i) = toc(t_start3);
    u23sJ = [boundry(length(t23sJ), a, t23sJ); u23sJ'];
    ODE23sJ_iter(i) = length(t23sJ) - 1;

end
Tbl = table(N, ODE23_iter, ODE23s_iter, ODE23sJ_iter, ODE23_time, ODE23s_time, ODE23sJ_time);
disp(Tbl)

%--------------------------------------
%functions 

function unext = expliciteuler(tn, un, A, bn, a, d, deltat, deltax)

if tn <= a
    bn(1) = d * sin(pi * tn / a) / deltax^2;
else 
    bn(1) = 0; 
end 

unext = un + deltat * (A * un + bn);

end 

function dydt = odefcn(t, u, A, b, a, d, deltax)

if t <= a
    b(1) = d * sin(pi * t / a) / deltax^2;
else 
    b(1) = 0; 
end 

dydt = A * u + b;

end 

function u_0_t = boundry(Nt, a, t)
u_0_t = zeros(1, Nt);
for j = 1:Nt
    if t(j) <= a
        u_0_t(j) = sin(pi * t(j) / a);
    else 
        u_0_t(j) = 0;
    end 
end 
end 

function [A, b, u, g] = Build_system(Nx, Nt, deltax, d, a)
A = zeros(Nx, Nx);
b = zeros(Nx, Nt);

for j = 1:Nx

    if j <= Nx - 1
        A(j, j) = -2 * d / deltax^2;
        A(j, j + 1) = d / deltax^2;
    end

    if j == Nx 
        A(j, j) = -2 * d / deltax^2;
        A(j, j - 1) = 2 * d / deltax^2;
    end

    if j >= 2 && j < Nx
        A(j, j - 1) = d / deltax^2;
    end

    u = zeros(Nx, Nt);
    g = zeros(Nx, 1);
end
A = sparse(A);
end
