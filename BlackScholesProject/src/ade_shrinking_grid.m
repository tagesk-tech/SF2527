function out = ade_shrinking_grid(cfg)
%ADE_SHRINKING_GRID Two-sweep Saul'yev/ADE Black--Scholes solver.
% Implements Wang et al. (2025), equations (8)--(10), while storing only
% the current time level. One rightmost point is discarded per time step.

cfg = defaults(cfg);
[x, u] = initial_grid_and_payoff(cfg);
initial_points = numel(x);
target_index = find(x <= cfg.S0, 1, 'last');
if initial_points - cfg.Nt < target_index + 1
    error('The shrinking grid loses the target stencil before final time.')
end

h = diff(x);
dt = cfg.T / cfg.Nt;
phi = zeros(size(u));
psi = zeros(size(u));
m = initial_points;
started = tic;

for n = 1:cfg.Nt
    mnew = m - 1;
    tau_new = n * dt;
    left = left_boundary(tau_new, cfg);
    phi(1) = left;

    for i = 2:mnew
        hL = h(i - 1);
        hR = h(i);
        hs = hL + hR;
        sx = cfg.sigma * x(i);
        fac_forward = 1 / dt + (sx^2 - cfg.r * x(i) * hR) / (hL * hs) ...
            - 0.5 * cfg.reaction_rate;
        numerator = u(i) / dt ...
            + sx^2 * ((u(i + 1) - u(i)) / (hR * hs) + phi(i - 1) / (hL * hs)) ...
            + cfg.r * x(i) * (hL * (u(i + 1) - u(i)) / (hR * hs) ...
            - hR * phi(i - 1) / (hL * hs)) ...
            + 0.5 * cfg.reaction_rate * u(i);
        phi(i) = numerator / fac_forward;
    end

    psi(mnew) = phi(mnew);
    for i = mnew - 1:-1:2
        hL = h(i - 1);
        hR = h(i);
        hs = hL + hR;
        sx = cfg.sigma * x(i);
        fac_backward = 1 / dt + (sx^2 + cfg.r * x(i) * hL) / (hR * hs) ...
            - 0.5 * cfg.reaction_rate;
        numerator = u(i) / dt ...
            + sx^2 * (psi(i + 1) / (hR * hs) - (u(i) - u(i - 1)) / (hL * hs)) ...
            + cfg.r * x(i) * (hL * psi(i + 1) / (hR * hs) ...
            + hR * (u(i) - u(i - 1)) / (hL * hs)) ...
            + 0.5 * cfg.reaction_rate * u(i);
        psi(i) = numerator / fac_backward;
    end
    psi(1) = left;
    u(1:mnew) = 0.5 * (phi(1:mnew) + psi(1:mnew));
    u(1) = left;
    m = mnew;
end

out = finish_output(x(1:m), u(1:m), cfg, toc(started), initial_points);
end

function cfg = defaults(cfg)
required = {'kind', 'K', 'r', 'sigma', 'T', 'Nt', 'h', 'Smax', 'S0'};
for k = 1:numel(required)
    if ~isfield(cfg, required{k})
        error('Missing cfg.%s', required{k})
    end
end
if ~isfield(cfg, 'C'), cfg.C = 100; end
if ~isfield(cfg, 'reaction_rate'), cfg.reaction_rate = -cfg.r; end
end

function [x, u] = initial_grid_and_payoff(cfg)
if strcmp(cfg.kind, 'cash-call') || strcmp(cfg.kind, 'cash-put')
    x = [0, 0.5 * cfg.h:cfg.h:cfg.Smax - 0.5 * cfg.h, cfg.Smax];
else
    n = round(cfg.Smax / cfg.h);
    x = linspace(0, cfg.Smax, n + 1);
end
switch cfg.kind
    case 'call'
        u = max(x - cfg.K, 0);
    case 'put'
        u = max(cfg.K - x, 0);
    case 'cash-call'
        u = cfg.C * double(x > cfg.K);
    case 'cash-put'
        u = cfg.C * double(x < cfg.K);
    otherwise
        error('Unknown option kind: %s', cfg.kind)
end
end

function value = left_boundary(tau, cfg)
switch cfg.kind
    case {'call', 'cash-call'}
        value = 0;
    case 'put'
        value = cfg.K * exp(-cfg.r * tau);
    case 'cash-put'
        value = cfg.C * exp(-cfg.r * tau);
end
end

function out = finish_output(x, u, cfg, runtime, initial_points)
[exact_profile, ~] = black_scholes_exact(x, cfg.T, cfg, cfg.kind);
out = struct();
out.x = x;
out.u = u;
out.exact = exact_profile;
out.target = target_value(x, u, cfg.S0);
out.target_exact = black_scholes_exact(cfg.S0, cfg.T, cfg, cfg.kind);
out.target_error = abs(out.target - out.target_exact);
out.max_error = max(abs(u - exact_profile));
out.runtime_seconds = runtime;
out.initial_points = initial_points;
out.final_points = numel(x);
out.final_right = x(end);
out.status = status_of(u);
end

function value = target_value(x, u, S0)
[distance, index] = min(abs(x - S0));
if distance <= 10 * eps(max(1, S0))
    value = u(index);
else
    value = interp1(x, u, S0, 'pchip');
end
end

function status = status_of(u)
if any(~isfinite(u))
    status = 'NONFINITE';
elseif max(abs(u)) > 1e8
    status = 'BLOWUP';
elseif min(u) < -1e-6
    status = 'NEGATIVE';
else
    status = 'OK';
end
end
