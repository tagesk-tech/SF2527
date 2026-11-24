function table_out = run_put_validation(root)
%RUN_PUT_VALIDATION Validate the group-requested European put extension.

Nt_values = [160; 320; 640; 1280];
put_price = zeros(size(Nt_values));
put_error = zeros(size(Nt_values));
parity_residual = zeros(size(Nt_values));
for k = 1:numel(Nt_values)
    put_cfg = project_config('put');
    put_cfg.Nt = Nt_values(k);
    put_out = ade_shrinking_grid(put_cfg);
    call_cfg = put_cfg;
    call_cfg.kind = 'call';
    call_out = ade_shrinking_grid(call_cfg);
    put_price(k) = put_out.target;
    put_error(k) = put_out.target_error;
    parity_residual(k) = abs(call_out.target - put_out.target ...
        - (put_cfg.S0 - put_cfg.K * exp(-put_cfg.r * put_cfg.T)));
end
observed_order = [NaN; log2(put_error(1:end-1) ./ put_error(2:end))];
table_out = table(Nt_values, put_price, put_error, observed_order, parity_residual);
writetable(table_out, fullfile(root, 'results', 'put_validation.csv'));
end

