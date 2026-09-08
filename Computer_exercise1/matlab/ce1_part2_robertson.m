%% SF2527 Computer Exercise 1 -- Part 2(a)
clear;
clc;
close all;

%% Robertson problem
r1 = 5.0e-2;
r2 = 1.2e4;
r3 = 4.0e7;

x0 = [1; 0; 0];
T = 10;

f = @(t, x) [-r1*x(1) + r2*x(2)*x(3);
              r1*x(1) - r2*x(2)*x(3) - r3*x(2)^2;
              r3*x(2)^2];

%% Runge-Kutta solution
h = 7.0e-4;
N = ceil(T/h);
h = T/N;

t = linspace(0, T, N + 1);
x = zeros(3, N + 1);
x(:, 1) = x0;

for n = 1:N
    x(:, n + 1) = rk3_step(f, t(n), x(:, n), h);
end

%% Concentrations on a linear scale
figure;
plot(t, x(1, :), 'LineWidth', 1.3);
hold on;
plot(t, x(2, :), 'LineWidth', 1.3);
plot(t, x(3, :), 'LineWidth', 1.3);
grid on;
xlabel('Time t');
ylabel('Concentration');
title('Robertson problem: explicit Runge-Kutta');
legend('x_A', 'x_B', 'x_C');

%% Concentrations on a logarithmic scale
figure;
semilogy(t, x(1, :), 'LineWidth', 1.3);
hold on;
semilogy(t, x(2, :), 'LineWidth', 1.3);
semilogy(t, x(3, :), 'LineWidth', 1.3);
grid on;
xlabel('Time t');
ylabel('Concentration');
title('Robertson problem: logarithmic scale');
legend('x_A', 'x_B', 'x_C');

fprintf('Stable RK step size: %.8f\n', h);



%% Part 2(b) Jacobian eigenvalues
J = @(x) [-r1,  r2*x(3), r2*x(2);
           r1, -r2*x(3) - 2*r3*x(2), -r2*x(2);
            0,  2*r3*x(2), 0];

lambda = zeros(3, N + 1);

for n = 1:N + 1
    lambda_n = eig(J(x(:, n)));
    [~, order] = sort(abs(lambda_n));
    lambda(:, n) = lambda_n(order);
end

figure;
semilogy(t, abs(lambda(2, :)), 'LineWidth', 1.3);
hold on;
semilogy(t, abs(lambda(3, :)), 'LineWidth', 1.3);
grid on;
xlabel('Time t');
ylabel('Eigenvalue magnitude');
title('Nonzero eigenvalues of the Jacobian');
legend('|\lambda_2|', '|\lambda_3|', 'Location', 'best');

R = @(z) 1 + z + z.^2/2 + z.^3/6;
z0 = fzero(@(z) abs(R(z)) - 1, [-3, -2]);
h_limit = abs(z0)/max(abs(lambda(3, :)));

fprintf('Largest Jacobian eigenvalue magnitude: %.6f\n', ...
        max(abs(lambda(3, :))));
fprintf('Predicted RK stability limit: %.8f\n', h_limit);


%% Part 2(c) Runge-Kutta solution on [0,1000]
T_long = 1000;
h_long = 2.8e-4;
N_long = ceil(T_long/h_long);
h_long = T_long/N_long;
x_long = x0;

tic;

for n = 1:N_long
    tn = (n - 1)*h_long;
    x_long = rk3_step(f, tn, x_long, h_long);
end

elapsed_time = toc;

fprintf('\nLong RK solution:\n');
fprintf('h = %.8f\n', h_long);
fprintf('x_A(1000) = %.12f\n', x_long(1));
fprintf('x_B(1000) = %.12e\n', x_long(2));
fprintf('x_C(1000) = %.12f\n', x_long(3));
fprintf('Computational time = %.4f seconds\n', elapsed_time);



%% Part 2(d) -- Implicit Euler
h_IE_tests = [10, 1, 0.1];

fprintf('\nImplicit Euler stability tests:\n');

for j = 1:length(h_IE_tests)
    [~, x_test] = impeuler(T, x0, h_IE_tests(j));
    fprintf('h = %.1f, x(10) = [%.6f, %.3e, %.6f]\n', ...
            h_IE_tests(j), x_test(end, 1), x_test(end, 2), x_test(end, 3));
end

h_IE = 0.1;
[t_IE, x_IE] = impeuler(T, x0, h_IE);

figure;
semilogy(t, x(1, :), '-', t_IE, x_IE(:, 1), '--', 'LineWidth', 1.3);
hold on;
semilogy(t, x(2, :), '-', t_IE, x_IE(:, 2), '--', 'LineWidth', 1.3);
semilogy(t, x(3, :), '-', t_IE, x_IE(:, 3), '--', 'LineWidth', 1.3);
grid on;
xlabel('Time t');
ylabel('Concentration');
title('Runge-Kutta and Implicit Euler');
legend('RK x_A', 'IE x_A', 'RK x_B', 'IE x_B', ...
       'RK x_C', 'IE x_C', 'Location', 'best');

[t_IE_long, x_IE_long] = impeuler(T_long, x0, h_IE);

figure;
semilogy(t_IE_long, x_IE_long(:, 1), 'LineWidth', 1.3);
hold on;
semilogy(t_IE_long, x_IE_long(:, 2), 'LineWidth', 1.3);
semilogy(t_IE_long, x_IE_long(:, 3), 'LineWidth', 1.3);
grid on;
xlabel('Time t');
ylabel('Concentration');
title('Implicit Euler solution on [0,1000]');
legend('x_A', 'x_B', 'x_C');




%% Part 2(e) Accuracy and efficiency
x_reference = [0.293414227164; 0.000001716342048; 0.706584056494];
h_IE_values = [1, 0.1, 0.01, 0.001];
IE_errors = zeros(length(h_IE_values), 1);
IE_times = zeros(length(h_IE_values), 1);

for j = 1:length(h_IE_values)
    tic;
    x_final = impeuler_nosave(T_long, x0, h_IE_values(j));
    IE_times(j) = toc;
    IE_errors(j) = norm(x_final - x_reference, 2);
end

method = [{'RK'}; repmat({'IE'}, length(h_IE_values), 1)];
step_size = [h_long; h_IE_values'];
endpoint_error = [norm(x_long - x_reference, 2); IE_errors];
computational_time = [elapsed_time; IE_times];

comparison = table(method, step_size, endpoint_error, computational_time);
disp(comparison);

figure;
loglog(IE_errors, IE_times, 'o-', 'LineWidth', 1.3);
hold on;
loglog(endpoint_error(1), computational_time(1), 's', 'MarkerSize', 8, ...
       'LineWidth', 1.3);
grid on;
xlabel('Error at t = 1000');
ylabel('Computational time (s)');
title('Efficiency of RK and Implicit Euler');
legend('Implicit Euler', 'Runge--Kutta', 'Location', 'best');
