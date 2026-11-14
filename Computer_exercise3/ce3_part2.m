clc; clearvars; close all
%% Part 2 parameters

Fr = 0.35;
alpha = 1/Fr;
L0 = -0.4;
L1 = 0.7;
T = 0.15;
Nx = 400;
requestedCFL = 0.9;

%% Grid, system matrix, and characteristic splitting
x = linspace(L0, L1, Nx)';
h = x(2) - x(1);

A = [1, alpha; alpha, 1];
[S, L] = eig(A);
Lp = L.*(L > 0);
Lm = L.*(L < 0);
Ap = S*Lp*S';
Am = S*Lm*S';

maximum_speed = max(abs(diag(L)));
[t, CFL] = make_time_grid(T, h, requestedCFL, maximum_speed);
k = t(2) - t(1);

f = @(x,t) (sin(40*pi*t + pi/6) > 0.5).*(abs(x) < 1/20).*sin(20*pi*x);

%% Part 2(b): analytical solution at the final time
[u_exact, v_exact] = analyticalsolution(x, T, alpha, f);

%% Part 2(c): numerical solutions on the same grid
[u_upwind, v_upwind] = upwindsystem(x, t, Ap, Am, f);
[u_LW, v_LW] = laxwendroffsystem(x, t, A, f);


%% Required final-time plots
figure(1)
plot(x, [u_exact, u_upwind(:,end), u_LW(:,end)], 'LineWidth', 1.2)
xlabel('x'), ylabel('u(x,0.15)'), title('Water level perturbation')
legend('Analytical', 'Upwind', 'Lax-Wendroff')
grid on

figure(2)
plot(x, [v_exact, v_upwind(:,end), v_LW(:,end)], 'LineWidth', 1.2)
xlabel('x'), ylabel('v(x,0.15)'), title('Water velocity perturbation')
legend('Analytical', 'Upwind', 'Lax-Wendroff')
grid on

%% 3D plots
u_exact_3D = zeros(length(x), length(t));
for n = 2:length(t)
    [u_exact_3D(:,n), ~] = analyticalsolution(x, t(n), alpha, f);
end

[X, Time] = meshgrid(x, t);

figure(3)
surf(X, Time, u_exact_3D')
xlabel('x'), ylabel('t'), zlabel('u(x,t)')
title('Analytical solution'), shading interp

figure(4)
surf(X, Time, u_upwind')
xlabel('x'), ylabel('t'), zlabel('u(x,t)')
title('Upwind solution'), shading interp

figure(5)
surf(X, Time, u_LW')
xlabel('x'), ylabel('t'), zlabel('u(x,t)')
title('Lax-Wendroff solution'), shading interp

%% Experiment: propagation speeds
% The slopes of the wave fronts show the two characteristic speeds.
figure(6)
imagesc(x, t, u_LW')
axis xy
colorbar
hold on
plot(0.05 + (1 + alpha)*t, t, 'k--', 'LineWidth', 1.2)
plot(-0.05 + (1 - alpha)*t, t, 'k--', 'LineWidth', 1.2)
hold off
xlabel('x'), ylabel('t'), title('Propagation of the two wave families')
legend('Right front: speed 1+alpha', 'Left front: speed 1-alpha')

%% Experiment: violate the CFL condition
% Compare the unstable solution growth with the stable run.
requestedCFL_unstable = 1.05;
[t_unstable, CFL_unstable] = make_time_grid( ...
    T, h, requestedCFL_unstable, maximum_speed);
assert(CFL_unstable > 1, 'The unstable experiment must violate the CFL limit.')
[u_unstable, v_unstable] = upwindsystem(x, t_unstable, Ap, Am, f);

size_stable = max([max(abs(u_upwind),[],1); max(abs(v_upwind),[],1)],[],1);
size_unstable = max([max(abs(u_unstable),[],1); max(abs(v_unstable),[],1)],[],1);

figure(7)
semilogy(t, size_stable, 'LineWidth', 1.2), hold on
semilogy(t_unstable, size_unstable, 'LineWidth', 1.2), hold off
xlabel('t'), ylabel('Maximum solution magnitude')
title('Stable and unstable upwind solutions')
legend(sprintf('CFL = %.3f', CFL), sprintf('CFL = %.3f', CFL_unstable), ...
    'Location', 'northwest')
grid on

%% Experiment: grid refinement at fixed CFL
% Compare the errors on successively finer grids.
Nx_vector = [50, 100, 200, 400];
h_vector = zeros(size(Nx_vector));
error_upwind = zeros(size(Nx_vector));
error_LW = zeros(size(Nx_vector));

for i = 1:length(Nx_vector)
    xi = linspace(L0, L1, Nx_vector(i))';
    hi = xi(2) - xi(1);
    [ti, ~] = make_time_grid(T, hi, CFL, maximum_speed);

    [ui_exact, ~] = analyticalsolution(xi, T, alpha, f);
    [ui_upwind, ~] = upwindsystem(xi, ti, Ap, Am, f);
    [ui_LW, ~] = laxwendroffsystem(xi, ti, A, f);

    h_vector(i) = hi;
    error_upwind(i) = max(abs(ui_upwind(:,end) - ui_exact));
    error_LW(i) = max(abs(ui_LW(:,end) - ui_exact));
end

figure(8)
loglog(h_vector, error_upwind, 'o-', 'LineWidth', 1.2), hold on
loglog(h_vector, error_LW, 'o-', 'LineWidth', 1.2), hold off
xlabel('Deltax'), ylabel('Maximum error in u(x,0.15)')
title('Grid refinement at fixed CFL number')
legend('Upwind', 'Lax-Wendroff', 'Location', 'northwest')
grid on

function [t, actualCFL] = make_time_grid(T, h, requestedCFL, maximum_speed)
% End exactly at T without exceeding the requested stable time step.
k_requested = requestedCFL*h/maximum_speed;
Nt = ceil(T/k_requested);
t = linspace(0, T, Nt + 1);
actualCFL = maximum_speed*(t(2) - t(1))/h;
end
