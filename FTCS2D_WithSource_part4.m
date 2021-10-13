function [uNt,E,tn] = FTCS2D_WithSource_part4(s,gl,gr,gb,gt,t0,u0,dt,Nt)

Lx = 10; Ly = 8;
Nx = 101; Ny = 80;
dx = Lx/(Nx-1);
dy = Ly/Ny;
X = ((0:Nx)-0.5)*dx;
Y = (0:Ny)*dy;
u2D = u0;
current_time = t0;
u2D(:,Nx+1) = u2D(:,Nx) - dx*gr(Y,current_time)';
u2D(:,1) = u2D(:,2) - dx*gl(Y,current_time)';
u2D(1,:) = gb(X,current_time);
u2D(Ny+1,:) = gt(X,current_time);


n = 1;
E(n) = 1;
tn(n) = current_time;
Err = 1;

    while ((n < Nt) && (Err > 10^-4 && Err < 10^3))
        u2DLast = u2D;
        u2D(2:Ny,2:Nx) = u2D(2:Ny,2:Nx)+...
            (dt/dx^2)*(u2D(2:Ny,3:Nx+1)-2*u2D(2:Ny,2:Nx)+u2D(2:Ny,1:Nx-1))...
           +(dt/dy^2)*(u2D(3:Ny+1,2:Nx)-2*u2D(2:Ny,2:Nx)+u2D(1:Ny-1,2:Nx))...
           +dt*s(u2D(2:Ny,2:Nx));

        current_time = n*dt;
        u2D(:,Nx+1) = u2D(:,Nx) + dx*gr(Y,current_time)';
        u2D(:,1) = u2D(:,2) - dx*gl(Y,current_time)';
        u2D(1,:) = gb(X,current_time);
        u2D(Ny+1,:) = gt(X,current_time);

        E(n) = max(max(abs((u2D(2:Ny,2:Nx)-u2DLast(2:Ny,2:Nx))/dt)));
        n = n + 1;
        tn(n) = current_time;
        Err = E(n-1);

    end
    tn = tn(1:end-1);
    uNt = u2D;
end