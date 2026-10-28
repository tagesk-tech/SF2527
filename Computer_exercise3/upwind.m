function u = upwind(x, t, a, g)
h = x(2) - x(1);
k = t(2) - t(1);
r = a*k/h; %CFL number, lambda = a*k/h
un = zeros(size(x));

for n = 1:length(t) - 1
    up = zeros(size(un));
    up(1) = g(t(n + 1));
    up(2:end) = un(2:end) - r*(un(2:end) - un(1:end-1));
    un = up;
end
u = un;
end
