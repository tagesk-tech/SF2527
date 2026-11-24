function tables = run_separated_convergence(root)
%RUN_SEPARATED_CONVERGENCE Refine time and space without moving final right edge.

final_right_target = 150;

Nt_values = [160; 320; 640; 1280];
target_error = zeros(size(Nt_values));
for k = 1:numel(Nt_values)
    cfg = project_config('call');
    cfg.Nt = Nt_values(k);
    cfg.Smax = final_right_target + cfg.Nt * cfg.h;
    out = ade_shrinking_grid(cfg);
    target_error(k) = out.target_error;
end
observed_order = [NaN; log2(target_error(1:end-1) ./ target_error(2:end))];
time_table = table(Nt_values, cfg.T ./ Nt_values, target_error, observed_order, ...
    'VariableNames', {'Nt', 'dt', 'target_error', 'observed_order'});
writetable(time_table, fullfile(root, 'results', 'time_convergence.csv'));

h_values = [0.4; 0.2; 0.1; 0.05];
Nt_fixed = 1280;
target_error = zeros(size(h_values));
max_error_near_target = zeros(size(h_values));
for k = 1:numel(h_values)
    cfg = project_config('call');
    cfg.Nt = Nt_fixed;
    cfg.h = h_values(k);
    cfg.Smax = final_right_target + cfg.Nt * cfg.h;
    out = ade_shrinking_grid(cfg);
    target_error(k) = out.target_error;
    mask = out.x >= 60 & out.x <= 140;
    max_error_near_target(k) = max(abs(out.u(mask) - out.exact(mask)));
end
target_order = [NaN; log2(target_error(1:end-1) ./ target_error(2:end))];
profile_order = [NaN; log2(max_error_near_target(1:end-1) ./ max_error_near_target(2:end))];
coupled_solver_table = table(h_values, target_error, target_order, ...
    max_error_near_target, profile_order);
writetable(coupled_solver_table, fullfile(root, 'results', 'space_coupled_solver_scan.csv'));

[operator_error, operator_order] = spatial_operator_study(h_values);
space_table = table(h_values, operator_error, operator_order);
writetable(space_table, fullfile(root, 'results', 'space_operator_convergence.csv'));

tables = struct('time', time_table, 'space_operator', space_table, ...
    'coupled_solver_scan', coupled_solver_table);
end

function [errors, orders] = spatial_operator_study(h_values)
% Check only the spatial Black--Scholes operator against exact u_tau.
p = project_config('call');
tau = 0.05;
errors = zeros(size(h_values));
Phi = @(z) 0.5 * erfc(-z ./ sqrt(2));
phi = @(z) exp(-0.5 * z.^2) ./ sqrt(2 * pi);
for k = 1:numel(h_values)
    h = h_values(k);
    x = 50:h:150;
    u = black_scholes_exact(x, tau, p, 'call');
    interior = 2:numel(x)-1;
    S = x(interior);
    ux = (u(interior + 1) - u(interior - 1)) / (2 * h);
    uxx = (u(interior + 1) - 2 * u(interior) + u(interior - 1)) / h^2;
    Lh = 0.5 * p.sigma^2 .* S.^2 .* uxx + p.r .* S .* ux - p.r .* u(interior);
    d1 = (log(S / p.K) + (p.r + 0.5 * p.sigma^2) * tau) / (p.sigma * sqrt(tau));
    d2 = d1 - p.sigma * sqrt(tau);
    u_tau = S .* phi(d1) * p.sigma / (2 * sqrt(tau)) ...
        + p.r * p.K * exp(-p.r * tau) .* Phi(d2);
    mask = S >= 60 & S <= 140;
    errors(k) = max(abs(Lh(mask) - u_tau(mask)));
end
orders = [NaN; log2(errors(1:end-1) ./ errors(2:end))];
end
