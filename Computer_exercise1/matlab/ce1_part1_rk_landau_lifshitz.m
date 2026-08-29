%% SF2527 Computer Exercise 1 -- Part 1(a)
clear;
clc;
close all;

%% Parameters and matrix A
alpha = 0.07;
a = (1/4) * [1; sqrt(11); 2];
m0 = [0; 0; 1];

C = [    0, -a(3),  a(2);
      a(3),     0, -a(1);
     -a(2),  a(1),     0];

A = C + alpha * C^2;
f = @(t, m) A * m;

%% Runge--Kutta method
T = 50;
N = 5000;
h = T/N;

t = linspace(0, T, N + 1);
m = zeros(3, N + 1);
m(:, 1) = m0;

for n = 1:N
    m(:, n + 1) = rk3_step(f, t(n), m(:, n), h);
end

%% Components of m(t)
figure;
plot(t, m(1, :), 'LineWidth', 1.3);
hold on;
plot(t, m(2, :), 'LineWidth', 1.3);
plot(t, m(3, :), 'LineWidth', 1.3);
grid on;
xlabel('Time t');
ylabel('Magnetization component');
title('Magnetization components');
legend('m_1', 'm_2', 'm_3');

%% Trajectory of m(t)
figure;
plot3(m(1, :), m(2, :), m(3, :), 'LineWidth', 1.4);
hold on;
plot3([0, a(1)], [0, a(2)], [0, a(3)], '--', 'LineWidth', 2);
grid on;
axis equal;
xlabel('m_1');
ylabel('m_2');
zlabel('m_3');
title('Magnetization trajectory');
legend('m(t)', 'a');

%% Order of accuracy
N_values = [50, 100, 200, 400, 800, 1600];
h_values = T ./ N_values;
mT = zeros(3, length(N_values));

for j = 1:length(N_values)
    N_test = N_values(j);
    h_test = T/N_test;
    u = m0;

    for n = 1:N_test
        tn = (n - 1)*h_test;
        u = rk3_step(f, tn, u, h_test);
    end

    mT(:, j) = u;
end

differences = zeros(1, length(N_values) - 1);

for j = 1:length(differences)
    differences(j) = norm(mT(:, j) - mT(:, j + 1), 2);
end

h_plot = h_values(1:end-1);
fit_points = 3:length(h_plot);
line_fit = polyfit(log(h_plot(fit_points)), log(differences(fit_points)), 1);
estimated_order = line_fit(1);

figure;
loglog(h_plot, differences, 'o-', 'LineWidth', 1.3);
hold on;
loglog(h_plot, exp(line_fit(2))*h_plot.^estimated_order, '--', 'LineWidth', 1.3);
grid on;
xlabel('Step size h');
ylabel('||m_N(T) - m_{2N}(T)||_2');
title('Runge-Kutta order of accuracy');
legend('Endpoint differences', sprintf('Slope = %.3f', estimated_order), ...
       'Location', 'northwest');

fprintf('Estimated order: %.6f\n', estimated_order);

%% Theoretical stability limit
eigenvalues = eig(A);
lambda = eigenvalues(imag(eigenvalues) > 0);
R = @(z) 1 + z + z.^2/2 + z.^3/6;
stability_equation = @(h) abs(R(h*lambda)) - 1;

h0 = fzero(stability_equation, [2, 2.5]);

fprintf('\nEigenvalues of A:\n');
disp(eigenvalues);
fprintf('Theoretical stability limit: %.6f\n', h0);

%% Empirical stability check
N_test = [ceil(T/h0), floor(T/h0)];
h_test = T ./ N_test;

figure;
hold on;

for j = 1:2
    t_test = linspace(0, T, N_test(j) + 1);
    u = zeros(3, N_test(j) + 1);
    u(:, 1) = m0;

    for n = 1:N_test(j)
        u(:, n + 1) = rk3_step(f, t_test(n), u(:, n), h_test(j));
    end

    plot(t_test, vecnorm(u, 2, 1), 'o-', 'LineWidth', 1.3);
end

grid on;
xlabel('Time t');
ylabel('||m(t)||_2');
title('Stable and unstable solutions');
legend(sprintf('Stable: h = %.4f', h_test(1)), ...
       sprintf('Unstable: h = %.4f', h_test(2)), ...
       'Location', 'northwest');

fprintf('Stable test:   h = %.6f, N = %d\n', h_test(1), N_test(1));
fprintf('Unstable test: h = %.6f, N = %d\n', h_test(2), N_test(2));
