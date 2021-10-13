function [uNt,E] = FTCS2D_WithSource(s,gl,gr,gb,gt,t0,u0,dt,Nt)


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

E = zeros(Nt,1);

    for n = 1:Nt
        u2DLast = u2D;
        u2D(2:Ny,2:Nx) = u2D(2:Ny,2:Nx)+...
            (dt/dx^2)*(u2D(2:Ny,3:Nx+1)-2*u2D(2:Ny,2:Nx)+u2D(2:Ny,1:Nx-1))...
           +(dt/dy^2)*(u2D(3:Ny+1,2:Nx)-2*u2D(2:Ny,2:Nx)+u2D(1:Ny-1,2:Nx))...
           +dt*s(u2D(2:Ny,2:Nx));

        current_time = t0 + n*dt;
        u2D(:,Nx+1) = u2D(:,Nx) + dx*gr(Y,current_time)';
        u2D(:,1) = u2D(:,2) - dx*gl(Y,current_time)';
        u2D(1,:) = gb(X,current_time);
        u2D(Ny+1,:) = gt(X,current_time);

        E(n) = max(max(abs((u2D(2:Ny,2:Nx)-u2DLast(2:Ny,2:Nx))/dt)));

    end
    uNt = u2D;
end