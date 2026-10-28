function [u, v] = analyticalsolution(x, t, alpha, f)
u = zeros(size(x));
v = zeros(size(x));
cp = 1 + alpha;
cm = 1 - alpha; %Characteristic speeds

for j = 1:length(x)
    fp = @(s) f(x(j) - cp*(t - s), s);
    fm = @(s) f(x(j) - cm*(t - s), s);

    Ip = integral(fp, 0, t, 'RelTol', 1e-10, 'AbsTol', 1e-12);
    Im = integral(fm, 0, t, 'RelTol', 1e-10, 'AbsTol', 1e-12);

    u(j) = (Ip - Im)/2;
    v(j) = (Ip + Im)/2;
end
end
