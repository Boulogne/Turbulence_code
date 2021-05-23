% Ex 1 flow in a channel with rotation
% Auteur : Boulogne Quentin
clear all;
close all
clc;
%% 1-Input Data
tic;
%Domain lengths
Lx=1;% L length along positive x axis [m]
Ly=0.5;% H hight along positive y axis [m]
domainLengths=[Lx Ly];% domain size [m]
%Mesh sizes
Nx=160; %Number of counter volumes on x axis
Ny=80;%Number of counter volumes on y axis
meshSizes=[Nx Ny];% Number of counter volumes on each axis
dx=Lx/(Nx+1);
dy=Ly/(Ny+1);
Npx=Nx+2;
Npy=Ny+2;
% Mesh definition
x=linspace(0,Lx,Npx);
y=linspace(Ly,0,Npy);
[X,Y] = meshgrid(x,y);
% Cylinder definition
solid = zeros(Npy,Npx);
Cylx=Npx/2;
Cyly=Npy/2;
rayon=0.1; %[m]
for i=1:Npy
for j=1:Npx
d=sqrt((dx*(j-Cylx))^2+(dy*(i-Cyly))^2);
if d<=rayon
solid(i,j)=1; % definition of the solid
end
end
end
% Fluid properties
Cp=1005;
R=287.058;
gamma=1.4;
% Initial flow conditions
Tin=300; % Temperature inflow [K]
Pin=1e5; % Pressure inflow [Pa]
Vin=20; % Initial Speed [m/s]
rhoin=Pin/R*Tin; % initial density [kg/m3]
% Computation data
error=1e-6; % error condition
max_number_iter=1e6;
%% 2-Initialisation of matrices
% Physics properties
Streamf=ones(Npy,Npx);
Temp=Tin*ones(Npy,Npx);
Pres=Pin*ones(Npy,Npx);
rho=rhoin*ones(Npy,Npx);
Vx=Vin*ones(Npy,Npx);
Vy=zeros(Npy,Npx);
V=Vx+Vy;
%% 3-Boundary conditions
Streamf(:,1)=Vin*y;
Streamf(1,:)=Vin*y(1);
Streamf(Npy,:)=0;
Streamfold=Streamf;
%% In the solid
for i=1:Npy
for j=1:Npx
if solid(i,j)==1
Streamf(i,j)=Vin*y(1)*3/4;% rotation condition 1/2 no rotation
Temp(i,j)=0;
Pres(i,j)=0;
rho(i,j)=0;
Vx(i,j)=0;
end
end
end
%% 4-Coefficients
rhoN = zeros(Npy,Npx);
rhoW = zeros(Npy,Npx);
rhoS = zeros(Npy,Npx);
rhoE = zeros(Npy,Npx);
aN=zeros(Npy,Npx);
aW=zeros(Npy,Npx);
aE=zeros(Npy,Npx);
aS=zeros(Npy,Npx);
aP=zeros(Npy,Npx);
bP=zeros(Npy,Npx);
errorini=1;
iteration=0;
while (errorini>error)
for i=2:Npy-1
for j=2:Npx-1
if solid(i,j)==0
rhoN(i,j)=dy/((dy/2)/(rhoin/rho(i,j))+(dy/2)/(rhoin/rho(i-1,j))); %harmonic mean
rhoS(i,j)=dy/((dy/2)/(rhoin/rho(i,j))+(dy/2)/(rhoin/rho(i+1,j)));
rhoE(i,j)=dx/((dx/2)/(rhoin/rho(i,j))+(dx/2)/(rhoin/rho(i,j+1))); %DENSITY
rhoW(i,j)=dx/((dx/2)/(rhoin/rho(i,j))+(dx/2)/(rhoin/rho(i,j-1)));
aN(i,j)=rhoN(i,j)*dx/dy;% discretization coefficient
aS(i,j)=rhoS(i,j)*dx/dy;
aE(i,j)=rhoE(i,j)*dy/dx;
aW(i,j)=rhoW(i,j)*dy/dx;
bP(i,j)=0; %in this case bp is equal to 0
aP(i,j)=aN(i,j)+aS(i,j)+aE(i,j)+aW(i,j);
Streamf(i,j) = (aE(i,j)*Streamf(i,j+1)+aW(i,j)*Streamf(i,j-1)+... % Stream fonction calculation
aN(i,j)*Streamf(i-1,j)+aS(i,j)*Streamf(i+1,j)+...
bP(i,j))/aP(i,j);
end
end
end
Streamf(:,Npx)=Streamf(:,Npx-1);% Outflow
errorini=max(max(abs(Streamf-Streamfold))); % convergence criteria
Streamfold=Streamf;% new stream function for next iteration
iteration=iteration+1;%Calculation of the number of iteration
end
%% 5-Computation of the other parameters
for i=2:Npy-1
for j=2:Npx-1
VxN=rhoN(i,j)*((Streamf(i-1,j)-Streamf(i,j))/dy);
VxS=rhoS(i,j)*((Streamf(i,j)-Streamf(i+1,j))/dy);
VyE=-rhoE(i,j)*((Streamf(i,j+1)-Streamf(i,j))/dx);
VyW=-rhoW(i,j)*((Streamf(i,j)-Streamf(i,j-1))/dx);
Vx(i,j)=(VxN+VxS)/2;
Vy(i,j)=(VyE+VyW)/2;
V(i,j)=sqrt(Vx(i,j)^2+Vy(i,j)^2); % Velocity [m/s]
Temp(i,j)=Tin+(Vin^2-V(i,j)^2)/(2*Cp); % Temperature [K]
Pres(i,j)=Pin*(Temp(i,j)/Tin)^(gamma/(gamma-1)); %Pressure calculation [Pa]
rho(i,j)=Pres(i,j)/(R*Temp(i,j)); % Density calculation [kg/m^3]
end
end
time_computation=toc;
%% 6-Post processing
C_p=zeros(Npy,Npx);
for i=1:Npy
for j=1:Npx
C_p(i,j) = 1-(V(i,j)/Vin)^2;
end
end
C_pmax=max(max(C_p));
figure(1)
contour(X,Y,C_p)
Mach=max(V)/345;
figure('Name','Streamfonctions','NumberTitle','off');
contour(X,Y,Streamf,30,'LineWidth',2);
title("Stream fonctions");
figure('Name','Flow velocity','NumberTitle','off');
quiver(X,Y,Vx,Vy);
title("Flow velocity");
figure;
pcolor(X,Y,V);
colorbar
shading interp;
axis equal
title('Velocity')
axis tight
grid on
figure;
pcolor(X,Y,rho);
colorbar
shading interp;
axis equal
title('Density')
axis tight
grid on
figure;
pcolor(X,Y,Pres);
colorbar
shading interp;
axis equal
title('Pressure')
axis tight
grid on
figure;
pcolor(X,Y,Temp);
colorbar
shading interp;
axis equal
title('Temperature')
axis tight
grid on