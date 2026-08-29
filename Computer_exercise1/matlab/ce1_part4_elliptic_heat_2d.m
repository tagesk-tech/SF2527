%% SF2527 Computer Exercise 1 -- Part 4
clear;
clc;
close all;

%% Parameters
Lx = 12;
Ly = 5;
Text = 25;
xpoint = 6;
ypoint = 2;

%% Part 4(a)
f = @(x, y) 2;
N = 60;
h = Lx/N;

[x, y, T] = solve_heat_2d(N, Lx, Ly, Text, f);

figure;
mesh(x, y, T);
xlabel('x');
ylabel('y');
zlabel('Temperature T');
title('Temperature for f = 2');

ix = round(xpoint/h) + 1;
iy = round(ypoint/h) + 1;
T_numerical = T(iy, ix);

fprintf('Numerical T(6,2) = %.12f\n', T_numerical);

%% Part 4(b) and Part 4(c)
c0 = Text;
c1 = 10;
c2 = -1;

T_exact = c0 + c1*ypoint + c2*ypoint^2;
error_constant = abs(T_numerical - T_exact);

fprintf('Exact T(6,2) = %.12f\n', T_exact);
fprintf('Error = %.3e\n', error_constant);

%% Part 4(d)
f = @(x, y) 100*exp(-0.5*(x - 4)^2 - 4*(y - 1)^2);
N_values = [60, 120, 240];
T_values = zeros(length(N_values), 1);

for j = 1:length(N_values)
    N = N_values(j);
    h = Lx/N;
    [x, y, T] = solve_heat_2d(N, Lx, Ly, Text, f);

    ix = round(xpoint/h) + 1;
    iy = round(ypoint/h) + 1;
    T_values(j) = T(iy, ix);

    fprintf('h = %.2f, T(6,2) = %.12f\n', h, T_values(j));

    if N == 120
        x_plot = x;
        y_plot = y;
        T_plot = T;
    end
end

observed_order = log2(abs(T_values(1) - T_values(2))/abs(T_values(2) - T_values(3)));
fprintf('Observed order: %.6f\n', observed_order);

figure;
mesh(x_plot, y_plot, T_plot);
xlabel('x');
ylabel('y');
zlabel('Temperature T');
title('Temperature from localized heat source');

figure;
imagesc(x_plot, y_plot, T_plot);
set(gca, 'YDir', 'normal');
axis equal tight;
cb = colorbar;
cb.Label.String = 'Temperature T';
xlabel('x');
ylabel('y');
title('Temperature from localized heat source');

figure;
contour(x_plot, y_plot, T_plot, 20);
axis equal tight;
cb = colorbar;
cb.Label.String = 'Temperature T';
xlabel('x');
ylabel('y');
title('Temperature from localized heat source');
