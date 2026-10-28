function [u, v] = laxwendroffsystem(x, t, A, f)
h = x(2) - x(1);
k = t(2) - t(1);
r = k/h; %lambda = k/h

u = zeros(length(x), length(t));
v = zeros(length(x), length(t));

for n = 1:length(t) - 1
    qn = [u(:,n)'; v(:,n)'];
    qp = zeros(size(qn)); %Solution at the next time level

    Fn = [zeros(1,length(x)); f(x,t(n))'];
    Fp = [zeros(1,length(x)); f(x,t(n+1))'];
    Ft = (Fn(:,2:end-1) + Fp(:,2:end-1))/2 - r*A*(Fn(:,3:end) - Fn(:,1:end-2))/4; %Corrected source

    qp(:,2:end-1) = qn(:,2:end-1) - r*A*(qn(:,3:end) - qn(:,1:end-2))/2 + r^2*A^2*(qn(:,3:end) - 2*qn(:,2:end-1) + qn(:,1:end-2))/2 + k*Ft;

    %Linear extrapolation to the two boundaries
    qp(:,1) = 2*qp(:,2) - qp(:,3);
    qp(:,end) = 2*qp(:,end-1) - qp(:,end-2);

    u(:,n+1) = qp(1,:)';
    v(:,n+1) = qp(2,:)';
end
end
