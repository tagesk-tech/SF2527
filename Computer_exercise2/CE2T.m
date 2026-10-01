clear; clc; close all


% initial values + function
Lx=12; Ly=5; Te=25; T=40; Nx=120; h=Lx/Nx; Ny=round(Ly/h);
dt=1/3; Nt=round(T/dt); dt=T/Nt; t=(0:Nt)*dt;
x=(0:Nx)*h; y=(0:Ny)*h; nx=Nx+1; n=nx*Ny;
f=@(x,y)100.*exp(-.5.*(x-4).^2-4.*(y-1).^2);


% creating matrix
A=spalloc(n,n,5*n); b=zeros(n,1);
for j=1:Ny
    for i=0:Nx
        r=i+1+(j-1)*nx;
        A(r,r)=-4/h^2;
        
        if i==0 %Neoman condition
            A(r,r+1)=2/h^2;
        elseif i==Nx %Neoman condition
            A(r,r-1)=2/h^2;
        
        else
            A(r,r-1)=1/h^2; 
            A(r,r+1)=1/h^2;
        end

        if j==1 % fixed lower boundrary
            b(r)=Te/h^2;
        else
            A(r,r-nx)=1/h^2;
        end
        if j==Ny % insulated upper boundrary
            A(r,r-nx)=A(r,r-nx)+1/h^2;
        else
            A(r,r+nx)=1/h^2;
        end
        b(r)=b(r)+f(i*h,j*h); % f(x_i, y_i)
    end
end

u=Te*ones(n,1); I=speye(n); 
L=I-dt*A/2; R=I+dt*A/2; % Lu^k+1 = Ru^k + dt*b
ts=[0 1 4 12 22 40]; ns=round(ts/dt); V=zeros(n,numel(ts)); V(:,1)=u;
rp=round(6/h)+1+(round(2/h)-1)*nx; % (x, y) = (6, 2)
z=zeros(Nt+1,1); z(1)=u(rp); 

for k=1:Nt
    u=L\(R*u+dt*b);
    z(k+1)=u(rp);
    for q=2:numel(ts)
        if k==ns(q)
            V(:,q)=u;
        end
    end
end

%% plot figures
for q=1:numel(ts)
    W=[Te*ones(nx,1),reshape(V(:,q),nx,Ny)];
    figure(q); mesh(x,y,W'); xlabel('x'); ylabel('y'); zlabel('u'); title(sprintf('t = %g',ts(q)))
end

figure(7); plot(t,z,'LineWidth',1.2); hold on; yline(47.224999298104,'--');
xlabel('t'); ylabel('u(6,2,t)'); legend('u(6,2,t)','steady value','Location','southeast'); grid on
fprintf('u(6,2,40) = %.10f\n',z(end))
