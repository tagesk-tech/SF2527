function [u, v] = upwindsystem(x, t, Ap, Am, f)
h = x(2) - x(1);
k = t(2) - t(1);
r = k/h; %lambda = k/h

u = zeros(length(x), length(t));
v = zeros(length(x), length(t));

for n = 1:length(t) - 1
    qn = [u(:,n)'; v(:,n)'];
    qp = zeros(size(qn)); %Solution at the next time level
    Fn = [zeros(1,length(x)); f(x,t(n))'];

    qp(:,2:end-1) = qn(:,2:end-1) - r*Ap*(qn(:,2:end-1) - qn(:,1:end-2)) - r*Am*(qn(:,3:end) - qn(:,2:end-1)) + k*Fn(:,2:end-1);

    %Linear extrapolation to the two boundaries
    qp(:,1) = 2*qp(:,2) - qp(:,3);
    qp(:,end) = 2*qp(:,end-1) - qp(:,end-2);

    u(:,n+1) = qp(1,:)';
    v(:,n+1) = qp(2,:)';
end
end
