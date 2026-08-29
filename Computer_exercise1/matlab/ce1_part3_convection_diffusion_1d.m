%% SF2527 Computer Exercise 1 -- Part 3(a)
clear;
clc;
close all;

%% Parameters
a = 0.1;
b = 0.4;
Q0 = 7000;
alpha0 = 50;
Tout = 25;
T0 = 100;

%% Part 3(a)
v = 1;
N_values = [10, 20, 40, 80];

figure;
hold on;

for j = 1:length(N_values)
    N = N_values(j);
    [z, T] = solve_convection_diffusion(N, v, a, b, Q0, alpha0, Tout, T0);
    plot(z, T, 'LineWidth', 1.3);
end

grid on;
xlabel('Position z');
ylabel('Temperature T');
title('Grid refinement for v = 1');
legend('N = 10', 'N = 20', 'N = 40', 'N = 80');

%% Temperature at z = 0.5
N_midpoint = [80, 160, 320];
T_midpoint = zeros(length(N_midpoint), 1);

for j = 1:length(N_midpoint)
    N = N_midpoint(j);
    [z, T] = solve_convection_diffusion(N, v, a, b, Q0, alpha0, Tout, T0);
    T_midpoint(j) = T(N/2 + 1);
    fprintf('N = %d, T(0.5) = %.12f\n', N, T_midpoint(j));
end

observed_order = log2(abs(T_midpoint(1) - T_midpoint(2))/abs(T_midpoint(2) - T_midpoint(3)));

fprintf('Observed order: %.6f\n', observed_order);

%% Part 3(b)
velocities = [1, 5, 15, 100];
N_velocity = [80, 100, 300, 2000];

figure;
hold on;

for j = 1:length(velocities)
    v = velocities(j);
    N = N_velocity(j);
    [z, T] = solve_convection_diffusion(N, v, a, b, Q0, alpha0, Tout, T0);
    plot(z, T, 'LineWidth', 1.3);
    fprintf('v = %g, N = %d, h = %.6f\n', v, N, 1/N);
end

grid on;
xlabel('Position z');
ylabel('Temperature T');
title('Temperature for different fluid velocities');
legend('v = 1', 'v = 5', 'v = 15', 'v = 100');
