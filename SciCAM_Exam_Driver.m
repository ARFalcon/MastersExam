clc; clear; close all;

alpha = 0;
s = @(u) alpha.*u;
f = @(x,y) 0.*x.*y;
gl = @(y,t) pi/10*(sinh(pi*(8-y)/5)/sinh(pi*8/5));
gr = @(y,t) pi/10*(sinh(pi*(8-y)/5)/sinh(pi*8/5));
gb = @(x,t) 1/2*sin(pi*x/5) + 1/2;
gt = @(x,t) -1/2*cos(pi*x/10) + 1/2;
t0 = 0;
dt = 0.002;
T = 160;
Nt = round(T/dt);
Lx = 10; Ly = 8;
Nx = 101; Ny = 80;
dx = Lx/(Nx-1);
dy = Ly/Ny;
X = ((0:Nx)-.5)*dx;
Y = (0:Ny)*dy;
[x2D,y2D] = meshgrid(X,Y);
u0 = f(x2D,y2D);
[u_sol,E] = FTCS2D_WithSource(s,gl,gr,gb,gt,t0,u0,dt,Nt);
tn = (t0:Nt-1)*dt;


%% Part 1 a) and b)

% a)
figure(1)
set(gca,'Fontsize',14)
hold on; box on; grid off;
set(gca,'Yscale','log')
semilogy(tn,E)
title('$E_n$ vs $t_n$ for $t=0$ to $t=160$','Interpreter','latex')
xlabel('$t_n$','Interpreter','latex')
ylabel('$E_n(t)$','Interpreter','latex')
saveas(gcf,'Part1a.png')


% b)
figure(2)
set(gca,'Fontsize',14)
hold on; box on; grid on;
surf(X,Y,u_sol,'edgecolor','none')
title('Numerical solution at t=160','Interpreter','latex')
xlabel('x','Interpreter','latex')
ylabel('y','Interpreter','latex')
zlabel('Numerical solution','Interpreter','latex')
view(3)
saveas(gcf,'Part1b.png')


%% Part 3 a), b) and c)
clc
Part3_Exact =@(x,y) 1/2*sin(pi.*x/5).*(sinh(pi.*(8-y)/5)/sinh(8*pi/5))-1/2*cos(pi.*x/10).*(sinh(pi.*y/10)/sinh(8*pi/10))+1/2;
Part3Exact = Part3_Exact(x2D,y2D);

% a)
figure(3)
set(gca,'Fontsize',14)
hold on; box on; grid on;
surf(X,Y,u_sol -Part3Exact,'edgecolor','none')
title('Difference between numerical solution and exact solution','Interpreter','latex')
xlabel('x','Interpreter','latex') 
ylabel('y','Interpreter','latex')
zlabel('Error of numerical solution','Interpreter','latex')
view(3)
saveas(gcf,'Part3a.png')


% b)
figure(4)
set(gca,'Fontsize',14)
hold on; box on; grid on;
plot(Y,u_sol(:,2)-Part3Exact(:,2),'-b','LineWidth',1.5)
plot(Y,u_sol(:,Nx)-Part3Exact(:,Nx),'-k','LineWidth',1.5)
xlabel('y')
ylabel('Difference')
title('Difference between $u^{N_t}_{k,j} - v(x_k,y_j)$','for different x-values','Interpreter','latex')
legend('$x=x_1$','$x=x_{N_x-1}$','Interpreter','latex')
saveas(gcf,'Part3b.png')

% c)

yt = 2/dy;
figure(5)
set(gca,'Fontsize',14)
hold on; box on; grid on;
plot(X,u_sol(yt,:)-Part3Exact(yt,:),'LineWidth',1.5)
xlabel('x','Interpreter','latex')
ylabel('Difference','Interpreter','latex')
title('Difference between $u^{N_t}_{i,k} - v(x_i,y_k)$',' at $y = 2$','Interpreter','latex')
saveas(gcf,'Part3c.png')


%% Part 4 Init



T = 800;
Nt = round(T/dt);
alpha = 0.1;
s = @(u) alpha.*u;
[u_sol_Part4_1,E_part4_1,tn_1] = FTCS2D_WithSource_part4(s,gl,gr,gb,gt,t0,u0,dt,Nt);

alpha = 0.14;
s = @(u) alpha.*u;
[u_sol_part4_2,E_part4_2,tn_2] = FTCS2D_WithSource_part4(s,gl,gr,gb,gt,t0,u0,dt,Nt);

alpha = 0.17;
s = @(u) alpha.*u;
[u_sol_part4_3,E_part4_3,tn_3] = FTCS2D_WithSource_part4(s,gl,gr,gb,gt,t0,u0,dt,Nt);

alpha = 0.2;
s = @(u) alpha.*u;
[u_sol_part4_4,E_part4_4,tn_4] = FTCS2D_WithSource_part4(s,gl,gr,gb,gt,t0,u0,dt,Nt);


%% Part 4 a), b) and c)

% a)

figure(6)
set(gca,'Fontsize',14)
hold on; box on; grid off;
set(gca,'Yscale','log')
semilogy(tn_1,E_part4_1,'LineWidth',1.5)
semilogy(tn_2,E_part4_2,'LineWidth',1.5)
semilogy(tn_3,E_part4_3,'LineWidth',1.5)
semilogy(tn_4,E_part4_4,'LineWidth',1.5)
xlabel('$t_n$','Interpreter','latex')
ylabel('$E_n$','Interpreter','latex')
title('$t_n$ vs $E_n$ for different alpha values','Interpreter','latex')
legend('$\alpha = 0.1$','$\alpha = 0.14$','$\alpha = 0.17$','$\alpha = 0.2$','location','E','Interpreter','latex')
saveas(gcf,'Part4a.png')

% b)

c_a1 = polyfit(tn_1(20:end),log(E_part4_1(20:end)),1);
c_a2 = polyfit(tn_2(20:end),log(E_part4_2(20:end)),1);
c_a3 = polyfit(tn_3(20:end),log(E_part4_3(20:end)),1);
c_a4 = polyfit(tn_4(20:end),log(E_part4_4(20:end)),1);

% c) 

alphavec = [0.1,0.14,0.17,0.2];
c1_vec = [c_a1(1),c_a2(1),c_a3(1),c_a4(1)];
figure(7)
set(gca,'Fontsize',14)
hold on; box on; grid on;
plot(alphavec,c1_vec)
xlabel('$\alpha$','Interpreter','latex')
ylabel('$c_1(\alpha)$','Interpreter','latex')
title('Plot of $c_1(\alpha)$ vs $\alpha$' ,'Interpreter','latex')
saveas(gcf,'Part4c.png')


%% Part 6 Init

eps = 10^-3;
s = @(u) 0.4*tanh(u);
f = @(x,y) eps*sin(pi*y/8);
gl = @(y,t) 0.*y;
gr = @(y,t) 0.*y;
gb = @(x,t) 0.*x;
gt = @(x,t) 0.*x;
t0 = 0;
dt = 0.002;
Lx = 10; Ly = 8;
Nx = 101; Ny = 80;
dx = Lx/(Nx-1);
dy = Ly/Ny;
X = ((0:Nx)-.5)*dx;
Y = (0:Ny)*dy;
[x2D,y2D] = meshgrid(X,Y);
u0 = f(x2D,y2D);
[u_sol_part6,E_part6, tn_part6] = FTCS2D_WithSource_part6(s,gl,gr,gb,gt,t0,u0,dt);


%% Part 6 a) and b)

% a)

figure(8)
set(gca,'Fontsize',14)
hold on; box on; grid off;
set(gca,'Yscale','log')
semilogy(tn_part6,E_part6,'LineWidth',1.5)
xlabel('$t_n$','Interpreter','latex')
ylabel('$E_n$','Interpreter','latex')
title('$t_n$ vs $E_n$ when $\epsilon = 10^{-3}$','Interpreter','latex')
saveas(gcf,'Part6a.png')

% b)

figure(9)
set(gca,'Fontsize',14)
surf(X,Y,u_sol_part6,'edgecolor','none')
xlabel('x','Interpreter','latex')
ylabel('y','Interpreter','latex')
zlabel('numerical solution','Interpreter','latex')
title('surface plot when $E(t_n)<10^{-8}$ for $\epsilon=10^{-3}$','Interpreter','latex')
saveas(gcf,'Part6b.png')


%% Part 7 Init


s = @(u) 0.1*tanh(u)-0.02;
f = @(x,y) 0.*x.*y;
gl = @(y,t) (pi*y/80).*cos(pi*t/10).^2;
gr = @(y,t) -(pi*y/80).*cos(pi*t/10).^2;
gb = @(x,t) 1/2*cos(pi*x/5).*cos(pi*t/5)+1/2;
gt = @(x,t) sin(pi*x/10).*cos(pi*t/10).^2;
t0 = 0;
dt = 0.002;
Lx = 10; Ly = 8;
Nx = 101; Ny = 80;
dx = Lx/(Nx-1);
dy = Ly/Ny;
X = ((0:Nx)-.5)*dx;
Y = (0:Ny)*dy;
[x2D,y2D] = meshgrid(X,Y);
u0 = f(x2D,y2D);
Tau = 10;
Np = round(Tau/dt);
k = 0;
t0 = k*Tau;
[u_sol_last,E_last] = FTCS2D_WithSource(s,gl,gr,gb,gt,t0,u0,dt,Np);
phi = 1;
while (k*Tau <= 400) && ((phi > 10^-8) && (phi < 10^3))
    u_sol_first = u_sol_last;
    k = k + 1;
    t0 = k*Tau;
    [u_sol_last,E_last] = FTCS2D_WithSource(s,gl,gr,gb,gt,t0,u_sol_first,dt,Np);
    Phi(k) = max(max(abs(u_sol_last(2:Ny,2:Nx)-u_sol_first(2:Ny,2:Nx))));
    phi = Phi(k);
end


%% Part 7 a), b) and c)

% a)

kt = (0:k-1)*Tau;
figure(10)
set(gca,'Fontsize',14)
hold on; box on; grid off;
set(gca,'Yscale','log')
semilogy(kt,Phi,'LineWidth',1.5)
xlabel('$k\tau$','Interpreter','latex')
ylabel('$\Phi_k$','Interpreter','latex')
title('Plot of $\Phi_k$ vs $k\tau$','Interpreter','latex')
saveas(gcf,'Part7a.png')

% b) 

figure(11)
set(gca,'Fontsize',14)
hold on; box on; grid on;
surf(X,Y,u_sol_last,'edgecolor','none')
xlabel('x','Interpreter','latex')
ylabel('y','Interpreter','latex')
zlabel('Numerical solution','Interpreter','latex')
title('Surface plot of final numerical solution','Interpreter','latex')
view(3)
saveas(gcf,'Part7b.png')

% c) 

yt_part7 = [4,6]/dy;
figure(12)
set(gca,'Fontsize',14)
hold on; box on; grid on;
plot(X,u_sol_last(yt_part7(1),:))
plot(X,u_sol_last(yt_part7(2),:))
xlabel('x','Interpreter','latex')
ylabel('numerical solution','Interpreter','latex')
title('Numerical solution at different y values','Interpreter','latex')
legend('$y=4$','$y=6$','Interpreter','latex')
saveas(gcf,'Part7c.png')
