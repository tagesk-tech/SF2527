function cfg = project_config(kind)
%PROJECT_CONFIG Parameters used in Wang et al. (2025).
cfg = struct('kind', kind, 'K', 100, 'r', 0.03, 'sigma', 0.3, ...
    'C', 100, 'T', 0.1, 'S0', 100, 'h', 0.05, ...
    'Smax', 400, 'Nt', 1280);
end

