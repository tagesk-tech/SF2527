function table_out = run_stability_stress(root)
%RUN_STABILITY_STRESS Compare ADE and Forward Euler, retaining failed cases.

Nt_values = [20; 40; 80; 160; 320; 500];
method = strings(2 * numel(Nt_values), 1);
Nt = zeros(size(method));
dt = zeros(size(method));
target_error = zeros(size(method));
max_abs_value = zeros(size(method));
min_value = zeros(size(method));
status = strings(size(method));
row = 0;

for k = 1:numel(Nt_values)
    cfg = project_config('call');
    cfg.h = 0.5;
    cfg.Smax = 400;
    cfg.Nt = Nt_values(k);
    outputs = {ade_shrinking_grid(cfg), forward_euler_shrinking_grid(cfg)};
    names = {'ADE', 'Forward Euler'};
    for j = 1:2
        row = row + 1;
        out = outputs{j};
        method(row) = names{j};
        Nt(row) = cfg.Nt;
        dt(row) = cfg.T / cfg.Nt;
        target_error(row) = out.target_error;
        max_abs_value(row) = max(abs(out.u));
        min_value(row) = min(out.u);
        status(row) = out.status;
    end
end

table_out = table(method, Nt, dt, target_error, max_abs_value, min_value, status);
writetable(table_out, fullfile(root, 'results', 'stability_stress.csv'));

fig = figure('Visible', 'off');
is_ade = method == "ADE";
semilogy(dt(is_ade), max_abs_value(is_ade), 'o-', 'LineWidth', 1.4); hold on
semilogy(dt(~is_ade), max_abs_value(~is_ade), 's-', 'LineWidth', 1.4);
xlabel('\Deltat'); ylabel('max |U|'); grid on
legend('ADE', 'Forward Euler', 'Location', 'best');
title('Stability stress test on the retained active domain');
exportgraphics(fig, fullfile(root, 'results', 'stability_stress.png'), 'Resolution', 180);
close(fig)
end

