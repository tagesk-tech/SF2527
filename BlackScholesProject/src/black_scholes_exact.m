function [price, greeks] = black_scholes_exact(S, tau, p, kind)
%BLACK_SCHOLES_EXACT Analytical Black--Scholes prices with no toolbox use.

arguments
    S double {mustBeNonnegative}
    tau (1,1) double {mustBeNonnegative}
    p struct
    kind (1,:) char
end

K = p.K;
r = p.r;
sigma = p.sigma;
C = getfield_default(p, 'C', 1);

if tau == 0
    switch kind
        case 'call'
            price = max(S - K, 0);
        case 'put'
            price = max(K - S, 0);
        case 'cash-call'
            price = C * double(S > K);
        case 'cash-put'
            price = C * double(S < K);
        otherwise
            error('Unknown option kind: %s', kind)
    end
    greeks = struct('delta', nan(size(S)), 'gamma', nan(size(S)));
    return
end

positive = S > 0;
d1 = -inf(size(S));
d2 = -inf(size(S));
d1(positive) = (log(S(positive) ./ K) + (r + 0.5 * sigma^2) * tau) ...
    ./ (sigma * sqrt(tau));
d2(positive) = d1(positive) - sigma * sqrt(tau);
Phi = @(z) 0.5 * erfc(-z ./ sqrt(2));
phi = @(z) exp(-0.5 * z.^2) ./ sqrt(2 * pi);
discount = exp(-r * tau);

switch kind
    case 'call'
        price = S .* Phi(d1) - K * discount .* Phi(d2);
        delta = Phi(d1);
        gamma = zeros(size(S));
        gamma(positive) = phi(d1(positive)) ./ (S(positive) * sigma * sqrt(tau));
    case 'put'
        price = K * discount .* Phi(-d2) - S .* Phi(-d1);
        delta = Phi(d1) - 1;
        gamma = zeros(size(S));
        gamma(positive) = phi(d1(positive)) ./ (S(positive) * sigma * sqrt(tau));
    case 'cash-call'
        price = C * discount .* Phi(d2);
        delta = nan(size(S));
        gamma = nan(size(S));
    case 'cash-put'
        price = C * discount .* Phi(-d2);
        delta = nan(size(S));
        gamma = nan(size(S));
    otherwise
        error('Unknown option kind: %s', kind)
end

greeks = struct('delta', delta, 'gamma', gamma);
end

function value = getfield_default(s, name, default)
if isfield(s, name)
    value = s.(name);
else
    value = default;
end
end

