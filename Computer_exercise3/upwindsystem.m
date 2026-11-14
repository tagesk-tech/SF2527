function [u, v] = upwindsystem(x, t, Ap, Am, f)
h = x(2) - x(1);
k = t(2) - t(1);
mu = k/h;

Nx = length(x);
Nt = length(t);
Q = zeros(2, Nx, Nt); %Q(:,j,n) = [u_j^n; v_j^n]

for n = 1:Nt - 1
    Qn = Q(:,:,n);
    Qnext = zeros(2, Nx);
    Fn = [zeros(1, Nx); f(x, t(n))'];

    backwardDifference = Qn(:,2:end-1) - Qn(:,1:end-2);
    forwardDifference = Qn(:,3:end) - Qn(:,2:end-1);

    Qnext(:,2:end-1) = Qn(:,2:end-1) ...
        - mu*Ap*backwardDifference ...
        - mu*Am*forwardDifference ...
        + k*Fn(:,2:end-1);

    Q(:,:,n+1) = extrapolate_system_boundaries(Qnext);
end

u = squeeze(Q(1,:,:));
v = squeeze(Q(2,:,:));
end
