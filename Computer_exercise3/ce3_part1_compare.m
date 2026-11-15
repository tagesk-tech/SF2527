clc, clear all, close all

T = 4;
D = 10;
a = 2;
tau = 2;
Nx = 100;
CFL = 0.9;

x = linspace(0, D, Nx)';
h = x(2) - x(1);
k = CFL*h/a;
Nt = ceil(T/k);
k = T/Nt;
t = (0:Nt)*k;

gsin = @(t) sin(2*pi*t/tau);
gsq = @(t) sign(sin(2*pi*t/tau));

u_exact = exactsolution(x, T, a, gsin);
u_upwind = upwind(x, t, a, gsin);
u_LF = laxfriedrichs(x, t, a, gsin);
u_LW = laxwendroff(x, t, a, gsin);

%% plots

figure(1)
plot(x, [u_exact, u_upwind, u_LF, u_LW], 'LineWidth', 1.2)
xlabel('x'), ylabel('u(x,4)'), title('Sine wave')
legend('Exact', 'Upwind', 'Lax-Friedrichs', 'Lax-Wendroff')
grid on

u_exact = exactsolution(x, T, a, gsq);
u_upwind = upwind(x, t, a, gsq);
u_LF = laxfriedrichs(x, t, a, gsq);
u_LW = laxwendroff(x, t, a, gsq);

figure(2)
plot(x, [u_exact, u_upwind, u_LF, u_LW], 'LineWidth', 1.2)
xlabel('x'), ylabel('u(x,4)'), title('Square wave')
legend('Exact', 'Upwind', 'Lax-Friedrichs', 'Lax-Wendroff')
grid on

%% d

g = gsin; %change to gsin for the sine wave
%Change upwind to laxfriedrichs or laxwendroff to test another method

%Different CFL numbers and fixed Nx
CFL_vector = [0.2, 0.4, 0.8, 0.9];

figure(3), hold on
for i = 1:length(CFL_vector)
    k = CFL_vector(i)*h/a;
    Nt = round(T/k);
    t = (0:Nt)*k;
    u = upwind(x, t, a, g);

    plot(x, u, 'LineWidth', 1.2)
end
hold off
xlabel('x'), ylabel('u(x,t)')
title('Upwind: different CFL numbers and fixed grid')
legend('CFL = 0.2', 'CFL = 0.4', 'CFL = 0.8', 'CFL = 0.9')
grid on

%Different Nx and fixed CFL number
Nx_vector = [50, 100, 200, 400];
CFL = 0.8;

figure(4), hold on
for i = 1:length(Nx_vector)
    x = linspace(0, D, Nx_vector(i))';
    h = x(2) - x(1);
    k = CFL*h/a;
    Nt = round(T/k);
    t = (0:Nt)*k;
    u = upwind(x, t, a, g);

    plot(x, u, 'LineWidth', 1.2)
end
hold off
xlabel('x'), ylabel('u(x,t)')
title('Upwind: different grids and CFL = 0.8')
legend('N = 50', 'N = 100', 'N = 200', 'N = 400')
grid on


%functions

function u = exactsolution(x, T, a, g)
tb = T - x/a;
u = zeros(size(x));
u(tb > 0) = g(tb(tb > 0));
end
