clc, clear all, close all
%Part 2

Fr = 0.35;
alpha = 1/Fr;
L0 = -0.4;
L1 = 0.7;
T = 0.15;
Nx = 400;
CFL = 0.9;
plot3D = true;

x = linspace(L0, L1, Nx)';
h = x(2) - x(1);

A = [1, alpha; alpha, 1];
[S, L] = eig(A);
Lp = L.*(L > 0);
Lm = L.*(L < 0);
Ap = S*Lp*inv(S);
Am = S*Lm*inv(S);

maximum_speed = max(abs(diag(L)));
k = CFL*h/maximum_speed;
Nt = ceil(T/k);
k = T/Nt;
t = (0:Nt)*k;

f = @(x,t) (sin(40*pi*t + pi/6) > 0.5).*(abs(x) < 1/20).*sin(20*pi*x);

%b: analytical solution at the final time
[u_exact, v_exact] = analyticalsolution(x, T, alpha, f);

%c: numerical solutions
[u_upwind, v_upwind] = upwindsystem(x, t, Ap, Am, f);
[u_LW, v_LW] = laxwendroffsystem(x, t, A, f);

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

if plot3D
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
end

%% Experiments

%The slopes of the wave fronts show the two propagation speeds
figure(6)
imagesc(x, t, u_LW')
axis xy, colorbar, hold on
plot(0.05 + (1 + alpha)*t, t, 'k--', 'LineWidth', 1.2)
plot(-0.05 + (1 - alpha)*t, t, 'k--', 'LineWidth', 1.2)
hold off
xlabel('x'), ylabel('t'), title('Propagation of the two wave families')
legend('Right front: speed 1+alpha', 'Left front: speed 1-alpha')

%Violate the CFL condition and compare the growth with the stable run
CFL_unstable = 1.05;
k_unstable = CFL_unstable*h/maximum_speed;
Nt_unstable = round(T/k_unstable);
k_unstable = T/Nt_unstable;
t_unstable = (0:Nt_unstable)*k_unstable;
[u_unstable, v_unstable] = upwindsystem(x, t_unstable, Ap, Am, f);

size_stable = max([max(abs(u_upwind),[],1); max(abs(v_upwind),[],1)],[],1);
size_unstable = max([max(abs(u_unstable),[],1); max(abs(v_unstable),[],1)],[],1);

figure(7)
semilogy(t, size_stable, 'LineWidth', 1.2), hold on
semilogy(t_unstable, size_unstable, 'LineWidth', 1.2), hold off
xlabel('t'), ylabel('Maximum solution magnitude')
title('Stable and unstable upwind solutions')
legend('CFL = 0.897', 'CFL = 1.05', 'Location', 'northwest')
grid on

%Compare the errors on successively finer grids
Nx_vector = [50, 100, 200, 400];
h_vector = zeros(size(Nx_vector));
error_upwind = zeros(size(Nx_vector));
error_LW = zeros(size(Nx_vector));

for i = 1:length(Nx_vector)
    xi = linspace(L0, L1, Nx_vector(i))';
    hi = xi(2) - xi(1);
    ki = CFL*hi/maximum_speed;
    Nti = ceil(T/ki);
    ki = T/Nti;
    ti = (0:Nti)*ki;

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
