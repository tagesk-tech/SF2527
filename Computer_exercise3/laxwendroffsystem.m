function [u, v] = laxwendroffsystem(x, t, A, f)
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
    Fnext = [zeros(1, Nx); f(x, t(n+1))'];

    centeredDifference = Qn(:,3:end) - Qn(:,1:end-2);
    secondDifference = Qn(:,3:end) ...
        - 2*Qn(:,2:end-1) + Qn(:,1:end-2);
    sourceAverage = (Fn(:,2:end-1) + Fnext(:,2:end-1))/2;
    sourceTransport = mu*A*(Fn(:,3:end) - Fn(:,1:end-2))/4;
    correctedSource = sourceAverage - sourceTransport;

    Qnext(:,2:end-1) = Qn(:,2:end-1) ...
        - mu*A*centeredDifference/2 ...
        + mu^2*A^2*secondDifference/2 ...
        + k*correctedSource;

    Q(:,:,n+1) = extrapolate_system_boundaries(Qnext);
end

u = squeeze(Q(1,:,:));
v = squeeze(Q(2,:,:));
end
