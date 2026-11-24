function tables = run_paper_reproduction(root)
%RUN_PAPER_REPRODUCTION Reproduce Tables 1 and 2 with corrected oracles.

Nt_values = [160; 320; 640; 1280];
kinds = {'call', 'cash-call'};
tables = struct();

for j = 1:numel(kinds)
    kind = kinds{j};
    price = zeros(size(Nt_values));
    corrected_exact = zeros(size(Nt_values));
    corrected_error = zeros(size(Nt_values));
    runtime_seconds = zeros(size(Nt_values));
    final_right = zeros(size(Nt_values));
    status = strings(size(Nt_values));
    for k = 1:numel(Nt_values)
        cfg = project_config(kind);
        cfg.Nt = Nt_values(k);
        out = ade_shrinking_grid(cfg);
        price(k) = out.target;
        corrected_exact(k) = out.target_exact;
        corrected_error(k) = out.target_error;
        runtime_seconds(k) = out.runtime_seconds;
        final_right(k) = out.final_right;
        status(k) = out.status;
        if Nt_values(k) == 1280
            make_profile_figure(out, kind, root);
        end
    end
    observed_order = [NaN; log2(corrected_error(1:end-1) ./ corrected_error(2:end))];
    table_out = table(Nt_values, price, corrected_exact, corrected_error, ...
        observed_order, runtime_seconds, final_right, status);
    tables.(strrep(kind, '-', '_')) = table_out;
    writetable(table_out, fullfile(root, 'results', ['paper_reproduction_' kind '.csv']));
end
end

function make_profile_figure(out, kind, root)
fig = figure('Visible', 'off');
mask = out.x <= 160;
plot(out.x(mask), out.exact(mask), 'k-', 'LineWidth', 1.8); hold on
plot(out.x(mask), out.u(mask), '--', 'LineWidth', 1.4);
xlabel('Asset price S'); ylabel('Option value'); grid on
legend('Analytical', 'ADE shrinking grid', 'Location', 'northwest');
title(strrep(kind, '-', ' '));
exportgraphics(fig, fullfile(root, 'results', ['profile_' kind '.png']), 'Resolution', 180);
close(fig)
end

