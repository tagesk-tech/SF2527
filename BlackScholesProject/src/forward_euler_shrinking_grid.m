function out = forward_euler_shrinking_grid(cfg)
%FORWARD_EULER_SHRINKING_GRID Independent explicit shrinking-grid baseline.

[x, u] = make_grid(cfg);
initial_points = numel(x);
target_index = find(x <= cfg.S0, 1, 'last');
if initial_points - cfg.Nt < target_index + 1
    error('The shrinking grid loses the target stencil before final time.')
end
h = diff(x);
dt = cfg.T / cfg.Nt;
m = initial_points;
started = tic;

for n = 1:cfg.Nt
    mnew = m - 1;
    unew = zeros(1, mnew);
    unew(1) = left_boundary(n * dt, cfg);
    for i = 2:mnew
        hL = h(i - 1);
        hR = h(i);
        hs = hL + hR;
        uxx = 2 * ((u(i + 1) - u(i)) / hR - (u(i) - u(i - 1)) / hL) / hs;
        ux = hL * (u(i + 1) - u(i)) / (hR * hs) ...
            + hR * (u(i) - u(i - 1)) / (hL * hs);
        Lu = 0.5 * cfg.sigma^2 * x(i)^2 * uxx + cfg.r * x(i) * ux - cfg.r * u(i);
        unew(i) = u(i) + dt * Lu;
    end
    u(1:mnew) = unew;
    m = mnew;
end

x = x(1:m);
u = u(1:m);
exact = black_scholes_exact(x, cfg.T, cfg, cfg.kind);
out = struct('x', x, 'u', u, 'exact', exact, ...
    'target', target_value(x, u, cfg.S0), ...
    'target_exact', black_scholes_exact(cfg.S0, cfg.T, cfg, cfg.kind), ...
    'runtime_seconds', toc(started), 'initial_points', initial_points, ...
    'final_points', m, 'final_right', x(end));
out.target_error = abs(out.target - out.target_exact);
out.max_error = max(abs(u - exact));
if any(~isfinite(u))
    out.status = 'NONFINITE';
elseif max(abs(u)) > 1e8
    out.status = 'BLOWUP';
elseif min(u) < -1e-6
    out.status = 'NEGATIVE';
else
    out.status = 'OK';
end

function value = target_value(x, u, S0)
[distance, index] = min(abs(x - S0));
if distance <= 10 * eps(max(1, S0))
    value = u(index);
else
    value = interp1(x, u, S0, 'pchip');
end
end
end

function [x, u] = make_grid(cfg)
n = round(cfg.Smax / cfg.h);
x = linspace(0, cfg.Smax, n + 1);
switch cfg.kind
    case 'call', u = max(x - cfg.K, 0);
    case 'put', u = max(cfg.K - x, 0);
    otherwise, error('Forward-Euler baseline currently supports call and put.')
end
end

function value = left_boundary(tau, cfg)
if strcmp(cfg.kind, 'put')
    value = cfg.K * exp(-cfg.r * tau);
else
    value = 0;
end
end
