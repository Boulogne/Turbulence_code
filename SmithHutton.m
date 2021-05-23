%Ex2 Smith-hutton problem
%Auteur Boulogne Quentin
clear all;
close all ;
clc;
%% INPUT DATA
clc;
clear;
%Domain lengths
Lx=2.0;% L length along positive x axis [m]
Ly=1.0;% H hight along positive y axis [m]
domainLengths=[Lx Ly];% domain size [m]
%Mesh sizes
Nx=199; %Number of counter volumes on x axis
Ny=199;%Number of counter volumes on y axis
meshSizes=[Nx Ny];% Number of counter volumes on each axis
scheme=1; %UDS = 1 CDS = 2 HDS = 3 EDS=4 PLDS=5
dx=Lx/(Nx+1);
dy=Ly/(Ny+1);
Npx=Nx+2;
Npy=Ny+2;
% Mesh definition
x=linspace(-1,1,Npx);
y=linspace(Ly,0,Npy);
[X,Y] = meshgrid(x,y);
an=zeros(Npy,Npx);
aw=zeros(Npy,Npx);
ae=zeros(Npy,Npx);
as=zeros(Npy,Npx);
ap=ones(Npy,Npx);
bp=zeros(Npy,Npx);
F_n=zeros(Npy,Npx);
F_s=zeros(Npy,Npx);
F_w=zeros(Npy,Npx);
F_e=zeros(Npy,Npx);
P_n=zeros(Npy,Npx);
P_s=zeros(Npy,Npx);
P_e=zeros(Npy,Npx);
P_w=zeros(Npy,Npx);
phi=zeros(Npy,Npx);
phinplus=zeros(Npy,Npx);
dt=0.001;
alpha=10;
S=0; % Source term
Coeff=1000000; % rho / diff
rho=1;% density [kg/m^3]
Diff_coeff=rho/Coeff;
Sx=dx;
Sy=dy;
Se=Sy;
Sw=Sy;
Ss=Sx;
Sn=Sx;
Vx=zeros(Npy,Npx);
Vy=zeros(Npy,Npx);
max_divergence=1.0e10;
max_time=10;
max_iteration=1000000;

tolerance_solver=1.0e-5;
error = 1;
for i=1:Npy
for j=1:Npx
De(i,j)=Diff_coeff*dy/dx;
Dw(i,j)=Diff_coeff*dy/dx;
Ds(i,j)=Diff_coeff*dx/dy;
Dn(i,j)=Diff_coeff*dx/dy;
Vx(i,j)=2*Y(i,j)*(1-(X(i,j))^2);
Vy(i,j)=-2*X(i,j)*(1-(Y(i,j))^2);
F_n(i,j)=rho*Vy(i,j)*Sn;
F_s(i,j)=rho*Vy(i,j)*Ss;
F_e(i,j)=rho*Vx(i,j)*Se;
F_w(i,j)=rho*Vx(i,j)*Sw;
P_n(i,j)=F_n(i,j)/Dn(i,j);
P_s(i,j)=F_s(i,j)/Ds(i,j);
P_e(i,j)=F_e(i,j)/De(i,j);
P_w(i,j)=F_w(i,j)/Dw(i,j);
if i > 1 && i < Npy && j > 1 && j < Npx && scheme==1 %UDS
ae(i,j)=De(i,j)*1+max(-F_e(i,j),0);
aw(i,j)=Dw(i,j)*1+max(F_w(i,j),0);
an(i,j)=Dn(i,j)*1+max(-F_n(i,j),0);
as(i,j)=Ds(i,j)*1+max(F_s(i,j),0);
ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+rho*dx*dy/dt;
end
if i > 1 && i < Npy && j > 1 && j < Npx && scheme==2 %CDS
ae(i,j)=De(i,j)*(1-0.5*abs(P_e(i,j)))+max(-F_e(i,j),0);
aw(i,j)=Dw(i,j)*(1-0.5*abs(P_w(i,j)))+max(F_w(i,j),0);
an(i,j)=Dn(i,j)*(1-0.5*abs(P_n(i,j)))+max(-F_n(i,j),0);
as(i,j)=Ds(i,j)*(1-0.5*abs(P_s(i,j)))+max(F_s(i,j),0);
ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+rho*dx*dy/dt;
end
if i > 1 && i < Npy && j > 1 && j < Npx && scheme==3 %HDS
ae(i,j)=De(i,j)*max(0,(1-0.5*abs(P_e(i,j))))+max(-F_e(i,j),0);
aw(i,j)=Dw(i,j)*max(0,(1-0.5*abs(P_w(i,j))))+max(F_w(i,j),0);
an(i,j)=Dn(i,j)*max(0,(1-0.5*abs(P_n(i,j))))+max(-F_n(i,j),0);
as(i,j)=Ds(i,j)*max(0,(1-0.5*abs(P_s(i,j))))+max(F_s(i,j),0);
ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+rho*dx*dy/dt;
end
if i > 1 && i < Npy && j > 1 && j < Npx && scheme==4 %EDS
ae(i,j)=De(i,j)*(abs(P_e(i,j))/(exp(abs(P_e(i,j)))-1))+max(-F_e(i,j),0);
aw(i,j)=Dw(i,j)*(abs(P_w(i,j))/(exp(abs(P_w(i,j)))-1))+max(F_w(i,j),0);
an(i,j)=Dn(i,j)*(abs(P_n(i,j))/(exp(abs(P_n(i,j)))-1))+max(-F_n(i,j),0);
as(i,j)=Ds(i,j)*(abs(P_s(i,j))/(exp(abs(P_s(i,j)))-1))+max(F_s(i,j),0);
ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+rho*dx*dy/dt;
end
if i > 1 && i < Npy && j > 1 && j < Npx && scheme==5 %PLDS
ae(i,j)=De(i,j)*max(0,(1-0.5*abs(P_e(i,j)))^5)+max(-F_e(i,j),0);
aw(i,j)=Dw(i,j)*max(0,(1-0.5*abs(P_w(i,j)))^5)+max(F_w(i,j),0);
an(i,j)=Dn(i,j)*max(0,(1-0.5*abs(P_n(i,j)))^5)+max(-F_n(i,j),0);
as(i,j)=Ds(i,j)*max(0,(1-0.5*abs(P_s(i,j)))^5)+max(F_s(i,j),0);
ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+rho*dx*dy/dt;
end
%max( P_n, P_s, P_e, P_w);
phi_old(i,j)=1-tanh(alpha);
% Inlet
if i==Npy && j<Npx/2+2
phi_old(i,j)=1+tanh(alpha*(2*X(i,j)+1));
end
%outlet
if i==Npy && j>Npx/2
phi_old(i,j)=phi_old(i-1,j);
end
end
end
phi=phi_old;
Nt=100000;
ite=0;
tend=0.01;
t=0;
while t<tend
for i=1:Npy
for j=1:Npx
bp(i,j)= S*dx*dy+rho*dx*dy*phi(i,j)/dt;
end
end
while error>tolerance_solver && ite<max_iteration
for i=1:Npy
for j=1:Npx
bp(i,j)= S*dx*dy+rho*dx*dy*phi(i,j)/dt;
if i==Npy && j<Npx/2+1
phinplus(i,j)=1+tanh(alpha*(2*X(i,j)+1));
end
if i==1
phinplus(i,j)=1-tanh(alpha);
end
if j==1 || j==Npx
phinplus(i,j)=1-tanh(alpha);
end
if i > 1 && i < Npy && j > 1 && j < Npx
phinplus(i,j)=(1/ap(i,j))*(ae(i,j)*phi(i,j+1)+aw(i,j)*phi(i,j-1)+an(i,j)*phi(i-1,j)+as(i,j)*phi(i+1,j)+bp(i,j));
end
if i==Npy && j> Npx/2
phinplus(i,j)=phinplus(i-1,j);
end
end
end
error=max(max(abs(phinplus-phi)));
phi=phinplus;
ite=ite+1;
end
t=t+dt;
phi=phinplus;
end
Stream=zeros(Npy,Npx);
for i=1:Npy
for j=1:Npx
Stream(i,j)=-(1-X(i,j)^2)*(1-Y(i,j)^2);
end
end
figure('Name','Streamfunction','NumberTitle','off');
contour(X,Y,Stream,10,'LineWidth',2);
title("Stream function");
Results=zeros(1,11);
Results(1)=phinplus(Npy,101);
Results(2)=phinplus(Npy,111);
Results(3)=phinplus(Npy,121);
Results(4)=phinplus(Npy,131);
Results(5)=phinplus(Npy,141);
Results(6)=phinplus(Npy,151);
Results(7)=phinplus(Npy,161);
Results(8)=phinplus(Npy,171);
Results(9)=phinplus(Npy,181);
Results(10)=phinplus(Npy,191);
Results(11)=phinplus(Npy,201);
Results;
% Analytical solution
x_ana=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
Coeff_10=[1.989 1.402 1.146 0.946 0.775 0.621 0.480 0.349 0.227 0.111 0.0];
Coeff_1000=[2.0 1.9990 1.9997 1.9850 1.8410 0.9510 0.1540 0.0010 0.0 0.0 0.0];
Coeff_1000000=[2.0 2.0 2.0 1.999 1.964 1.000 0.036 0.001 0.0 0.0 0.0];
Coeff_anal=[Coeff_10 Coeff_1000 Coeff_1000000];
% plot(x_ana,phinplus(Npy,1:11));
figure('Name','Analatical results results ','NumberTitle','off');
plot(x_ana,Coeff_10,'-o',x_ana,Coeff_1000,'-o',x_ana,Coeff_1000000,'-o','LineWidth',2)
legend('coeff = 10','coeff=1000','coeff=1e6')
title('Analytical Results')
figure('Name','Result 2D vizualisation','NumberTitle','off');
pcolor(X,Y,phinplus);
colorbar
shading interp;
axis equal
title('Velocity repartition')
axis tight
grid on
totalresult=phinplus(Npy,101:end);
if Coeff==10
error_result=max(abs(Results-Coeff_10)*100/Coeff_10);
figure('Name','Results Comparaison ','NumberTitle','off');
plot(x_ana,Coeff_10,'-o',x_ana,Results,'-o','LineWidth',2)
legend('coeff = 10','Numerical results')
figure('Name','Results for the mesh used ','NumberTitle','off');
plot(X(Npy,101:end),totalresult,'LineWidth',2)
legend('coeff = 10')
xlabel('x');
ylabel('v');
grid
end
if Coeff==1000
error_result=max(abs(Results-Coeff_1000)*100/Coeff_1000);
figure('Name','Results Comparaison ','NumberTitle','off');
plot(x_ana,Coeff_1000,'-o',x_ana,Results,'-o','LineWidth',2)
legend('coeff = 1000','Numerical results')
figure('Name','Results for the mesh used ','NumberTitle','off');
plot(X(Npy,101:end),totalresult,'LineWidth',2)
legend('rho/diff = 1000')
xlabel('x');
ylabel('v');
grid
end
if Coeff==1000000
error_result=max(abs(Results-Coeff_1000000)*100/Coeff_1000000);
figure('Name','Results Comparaison ','NumberTitle','off');
plot(x_ana,Coeff_1000000,'-o',x_ana,Results,'-o','LineWidth',2)
legend('coeff = 1000000','Numerical results')
figure('Name','Results for the mesh used ','NumberTitle','off');
plot(X(Npy,101:end),totalresult,'LineWidth',2)
xlabel('x');
ylabel('v');
legend('rho/diff = 1000000')
grid
end
Peclet=max(max(abs(P_n)));