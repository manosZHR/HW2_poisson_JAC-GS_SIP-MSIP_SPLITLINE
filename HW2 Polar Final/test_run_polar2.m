clear; clc
close all
Nr=50;
Nt=50;
itgmr=200000;
omega=0.4;
psi=0.99;
rorigin=1;
rlast=2;
torigin=0;
tlast=2*pi();
tol=1e-10;


[u_gs,res_gs,iter_gs,R_gs, T_gs] = GaussSeidelPolarNeumann(Nr,Nt,rorigin,rlast,torigin,tlast,itgmr,1,tol);
[u_jac,res_jac,iter_jac] = JacobiPolarNeumann(Nr,Nt,rorigin,rlast,torigin,tlast,itgmr,omega,tol);
[u_sip5n,res_sip5d,iter_sip5] = MSIP5PolarNeumann(Nr,Nt,rorigin,rlast,torigin,tlast,0,itgmr,tol);
[u_sip9n,res_sip9d,iter_sip9] = MSIP9PolarNeumann(Nr,Nt,rorigin,rlast,torigin,tlast,0,itgmr,tol);
[u_msip5n,res_msip5,iter_msip5] = MSIP5PolarNeumann(Nr,Nt,rorigin,rlast,torigin,tlast,psi,itgmr,tol);
[u_msip9n,res_msip9,iter_msip9] = MSIP9PolarNeumann(Nr,Nt,rorigin,rlast,torigin,tlast,psi,itgmr,tol);

[u_gsd,res_gsd,iter_gsd,R_gs, T_gs] = GaussSeidelPolarDirichlet(Nr,Nt,rorigin,rlast,torigin,tlast,itgmr,1,tol);
[u_jacd,res_jacd,iter_jacd] = JacobiPolarDirichlet(Nr,Nt,rorigin,rlast,torigin,tlast,itgmr,1,tol);
[u_sip5d,res_sip5d,iter_sip5d] = MSIP5PolarDirichlet(Nr,Nt,rorigin,rlast,torigin,tlast,0,itgmr,tol);
[u_sip9d,res_sip9d,iter_sip9d] = MSIP9PolarDirichlet(Nr,Nt,rorigin,rlast,torigin,tlast,0,itgmr,tol);
[u_msip5d,res_msip5d,iter_msip5d] = MSIP5PolarDirichlet(Nr,Nt,rorigin,rlast,torigin,tlast,psi,itgmr,tol);
[u_msip9d,res_msip9d,iter_msip9d] = MSIP9PolarDirichlet(Nr,Nt,rorigin,rlast,torigin,tlast,psi,itgmr,tol);



Rp = reshape(R_gs,Nr,Nt); %check if reshape is correct
Tp = reshape(T_gs,Nr,Nt);

figure(1)
up = reshape(u_msip9d,Nr,Nt);

hold on
box on
xlabel('x', 'FontSize', 20),ylabel('y', 'FontSize', 20),title('Surface plot of u(x,y) (Dririchlet bc)', 'FontSize', 20)
surf(Rp.*cos(Tp),Rp.*sin(Tp),up) %converting to cartesian coordinates
colorbar
hold off

figure(2)
up = reshape(u_msip9n,Nr,Nt);

hold on
box on
xlabel('x', 'FontSize', 20),ylabel('y', 'FontSize', 20),title('Surface plot of u(x,y)(Neumann bc)', 'FontSize', 20)
surf(Rp.*cos(Tp),Rp.*sin(Tp),up) %converting to cartesian coordinates
colorbar
hold off