root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'), fullfile(root, 'experiments'));

p = project_config('call');
call_exact = black_scholes_exact(100, 0.1, p, 'call');
cash_exact = black_scholes_exact(100, 0.1, p, 'cash-call');
assert(abs(call_exact - 3.929276) < 1e-6, 'Correct call oracle failed.')
assert(abs(cash_exact - 49.221357) < 1e-6, 'Cash-call oracle failed.')
assert(abs(call_exact - 3.9123) > 1e-2, 'Paper benchmark defect was not exposed.')

cfg = project_config('call');
cfg.Nt = 320;
out = ade_shrinking_grid(cfg);
assert(abs(out.target - 3.7501) < 7e-3, 'Paper table reproduction failed.')

cfg.Nt = 640;
fine = ade_shrinking_grid(cfg);
assert(fine.target_error < out.target_error, 'Refinement did not reduce target error.')

put_cfg = project_config('put');
put_cfg.Nt = 640;
put_out = ade_shrinking_grid(put_cfg);
parity = abs(fine.target - put_out.target ...
    - (put_cfg.S0 - put_cfg.K * exp(-put_cfg.r * put_cfg.T)));
assert(parity < 2e-2, 'Numerical put--call parity failed.')
assert(abs(put_out.u(1) - put_cfg.K * exp(-put_cfg.r * put_cfg.T)) < 1e-12, ...
    'Put left boundary failed.')

fault_cfg = project_config('call');
fault_cfg.Nt = 1280;
correct = ade_shrinking_grid(fault_cfg);
fault_cfg.reaction_rate = +fault_cfg.r; % planted wrong sign on the reaction term
faulty = ade_shrinking_grid(fault_cfg);
assert(faulty.target_error > correct.target_error && abs(faulty.target - correct.target) > 2e-2, ...
    'Analytical oracle did not reject the planted reaction-sign fault.')

fprintf('All tests passed. Planted reaction-sign fault was rejected.\n')
