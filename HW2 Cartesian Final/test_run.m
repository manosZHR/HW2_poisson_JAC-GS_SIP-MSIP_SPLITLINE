clear; close all; clc
Nx=50;
Ny=50;
itgmr=50000;
omega=1;
psi=0.8;
xmin=-2;
xmax=2;
ymin=-2;
ymax=2;
tol=1e-10;

[u_jac,res_jac,iter_jac] = Jacobi(Nx,Ny,xmin,xmax,ymin,ymax,itgmr,omega,tol);
[u_gs,res_gs,iter_gs] = GaussSeidel(Nx,Ny,xmin,xmax,ymin,ymax,itgmr,omega,tol);
[u_sip5,res_sip5,iter_sip5] = MSIP5(Nx,Ny,xmin,xmax,ymin,ymax,0,itgmr,tol);
[u_sip9,res_sip9,iter_sip9] = MSIP9(Nx,Ny,xmin,xmax,ymin,ymax,0,itgmr,tol);
[u_msip5,res_msip5,iter_msip5] = MSIP5(Nx,Ny,xmin,xmax,ymin,ymax,psi,itgmr,tol);
[u_msip9,res_msip9,iter_msip9] = MSIP9(Nx,Ny,xmin,xmax,ymin,ymax,psi,itgmr,tol);

dx = (xmax-xmin)/(Nx-1);
dy = (ymax-ymin)/(Nx-1);
x = xmin:dx:xmax;
y = ymin:dy:ymax;
[xcoord,ycoord] = meshgrid(x,y);
ycoord=flip(ycoord);

up = [];

for j=1:Ny

    for i = 1:Nx

        up(i,j) = u_msip9(j*Ny+1-i);

    end
    
end


f1 = figure(1);
hold on
box on
xlabel('x', 'FontSize', 20),ylabel('y', 'FontSize', 20),title('Contour plot of the surface u(x,y)', 'FontSize', 20)
contour(xcoord,ycoord,up)
hold off


