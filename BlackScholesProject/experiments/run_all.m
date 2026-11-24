root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'src'), fullfile(root, 'experiments'));
if ~exist(fullfile(root, 'results'), 'dir'), mkdir(fullfile(root, 'results')); end

fprintf('Running paper reproduction...\n')
paper_tables = run_paper_reproduction(root);
disp(paper_tables.call)
disp(paper_tables.cash_call)

fprintf('Running separated convergence studies...\n')
convergence_tables = run_separated_convergence(root);
disp(convergence_tables.time)
disp(convergence_tables.space_operator)
disp(convergence_tables.coupled_solver_scan)

fprintf('Running put validation...\n')
put_table = run_put_validation(root);
disp(put_table)

fprintf('Running stability stress test...\n')
stability_table = run_stability_stress(root);
disp(stability_table)

save(fullfile(root, 'results', 'analysis_workspace.mat'), 'paper_tables', ...
    'convergence_tables', 'put_table', 'stability_table');
fprintf('Analysis complete. Results are in %s\n', fullfile(root, 'results'))
