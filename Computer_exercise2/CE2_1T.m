clear; clc; close all


% initial values
T=2; d=.35; a=1.2; tau1=1.1;
Nx=50; dx=1/Nx;

dtmax=dx^2/(2*d);             % explicit Euler stability limit
dt=.99*dtmax;
[t,x,U,dt]=explicit_euler(T,d,a,Nx,dt);


% plots
figure(1); surf(t,x,U,'EdgeColor','none')
xlabel('\tau'); ylabel('x'); zlabel('u'); title('u(x,\tau)')

u1=interp1(t,U',tau1)';
figure(2); plot(x,u1,'LineWidth',1.2); grid on
xlabel('x'); ylabel('u'); title(sprintf('u(x,%.1f)',tau1))

figure(3); plot(t,U(1,:),'LineWidth',1.2); hold on
plot(t,U(end,:),'LineWidth',1.2); grid on
xlabel('\tau'); ylabel('u'); legend('u(0,\tau)','u(1,\tau)','Location','best')


% stability comparison
[tu,xu,Uu,dtu]=explicit_euler(T,d,a,Nx,1.01*dtmax);
fprintf('Stable:   d*dt/dx^2 = %.6f\n',d*dt/dx^2)
fprintf('Unstable: d*dt/dx^2 = %.6f\n',d*dtu/dx^2)

figure(4); surf(tu,xu,Uu,'EdgeColor','none')
xlabel('\tau'); ylabel('x'); zlabel('u'); title('Unstable solution')


% part b)
N=[100;200;400];
tspan=[0 T];
options1=odeset('RelTol',1e-5);

ODE23_time=zeros(length(N),1);
ODE23s_time=zeros(length(N),1);
ODE23sJ_time=zeros(length(N),1);

ODE23_iter=zeros(length(N),1);
ODE23s_iter=zeros(length(N),1);
ODE23sJ_iter=zeros(length(N),1);

for i=1:length(N)
    dx=1/N(i); k=d/dx^2;
    A=-2*eye(N(i))+diag(ones(N(i)-1,1),1)+diag(ones(N(i)-1,1),-1);
    A(end,end-1)=2;
    A=k*sparse(A);

    u0=zeros(N(i),1);
    rhs=@(tau,u) odefcn(tau,u,A,k,a);

    tic
    sol23=ode23(rhs,tspan,u0,options1);
    ODE23_time(i)=toc;
    ODE23_iter(i)=length(sol23.x)-1;

    tic
    sol23s=ode23s(rhs,tspan,u0,options1);
    ODE23s_time(i)=toc;
    ODE23s_iter(i)=length(sol23s.x)-1;

    options2=odeset(options1,'Jacobian',A);
    tic
    sol23sJ=ode23s(rhs,tspan,u0,options2);
    ODE23sJ_time(i)=toc;
    ODE23sJ_iter(i)=length(sol23sJ.x)-1;
end

Tbl=table(N,ODE23_iter,ODE23s_iter,ODE23sJ_iter,...
    ODE23_time,ODE23s_time,ODE23sJ_time);
disp(Tbl)


% function
function [t,x,U,dt]=explicit_euler(T,d,a,Nx,dt)
dx=1/Nx;
Nt=ceil(T/dt); dt=T/Nt;
x=(0:Nx)*dx; t=(0:Nt)*dt;

k=d/dx^2;
A=-2*eye(Nx)+diag(ones(Nx-1,1),1)+diag(ones(Nx-1,1),-1);
A(end,end-1)=2;
A=k*sparse(A);

g=zeros(1,Nt+1);
heated=t<=a;
g(heated)=sin(pi*t(heated)/a);

b=zeros(Nx,Nt+1);
b(1,:)=k*g;

u=zeros(Nx,Nt+1);
for n=1:Nt
    u(:,n+1)=u(:,n)+dt*(A*u(:,n)+b(:,n));
end
U=[g;u];
end

function dudt=odefcn(t,u,A,k,a)
g=0;
if t<=a
    g=sin(pi*t/a);
end
b=zeros(length(u),1);
b(1)=k*g;
dudt=A*u+b;
end
