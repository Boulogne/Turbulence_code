%eX4 bURGER EQUATION
%aUTEUR bOULOGNE qUENTIN
clc;
clear ;
close all;
%% Input data
Re = 40;
N = 20;
C1 = 0.01;
Ck = 0.05; %Kolmogorov constant
m=2;
LES = 0; %0 desactivation of LES (Turbulence) 1 Activation
dt=(C1*Re)/(N^2);
Toler_conv=1e-6;
%% Without turbulence
LES = 0;
Ck = 0.05;
E1=CalcultationofBurgerequation(Ck,LES,N,Re,dt,m,Toler_conv,0.001);
%% LES turbulence Ck=0.4223
LES = 1;
Ck = 0.4223;
E2=CalcultationofBurgerequation(Ck,LES,N,Re,dt,m,Toler_conv,0.001);
%% LES turbulence Ck=0.05
LES = 1;
Ck = 0.05;
E3=CalcultationofBurgerequation(Ck,LES,N,Re,dt,m,Toler_conv,0.001);
%% N=100
LES = 0;
Ck = 0.05;
E4=CalcultationofBurgerequation(Ck,LES,100,Re,dt,m,1e-5,0.001);
N=100;
for k=1:N
K(1,k)=k;
Slp_1(k)=1/(k^2);
end
figure(1);
loglog(K,E4,'-+',K,Slp_1);
grid;
xlabel('k');
ylabel('E_k');
legend('N=100 DNS','slope=-2');
%% Post treatement
N=20;
for k=1:N
K20(1,k)=k;
Slp_2(k)=1/(k^2);
end
figure(2);
loglog(K20,E1,'-+',K20,E2,'-+',K20,E3,'-+',K20,Slp_2);
grid;
xlabel('k');
ylabel('E_k');
legend('N=20 DNS','N=20 LES Ck=0.4223','N=20 LES Ck=0.05','slope=-2');


function E = CalcultationofBurgerequation(Ck,LES,N,Re,dt,m,Toler_conv,C1)
Re=40;
dt=(C1*Re)/(N^2);
u0 = zeros(N,1);
K = zeros(N,1);
for k = 1:N
K(k,1) = k;
u0(k) =1/abs(k);
end
t = 1;
u(:,t) =u0(:);
Diffusion=zeros(N,1);
Convection=zeros(N,1);
for k=1:N
if LES==0 % Without LES
Diffusion(k,1)=k^2*(1/Re)*u(k);
end
if LES==1% With LES
Vt_infinite=0.31*(5-m)/(m+1)*sqrt(3-m)*Ck^(-3/2);
Vt_star=1+34.5*exp(-3.03*N/k);
Ekn=u(N)*conj(u(N));
Vt=Vt_infinite*sqrt(Ekn/N)*Vt_star;
Diffusion(k,1)=k^2*(1/Re+Vt)*u(k);
end
for p=-N+k:N
q=k-p;
if q~=0 && p~=0
if q<0
uq=conj(u(-q));
else
uq=u(q);
end
if p<0
up=conj(u(-p));
else
up=u(p);
end
Convection(k)=Convection(k)+up*q*1i*uq;
end
end
end
Sum_ini(:,1)=Diffusion(K)+Convection(K);
t = t+1;
u(:,t) = u(:,t-1);
%% Solution
converging_condition = 0;
while converging_condition == 0
Diffusion=zeros(N,1);
Convection=zeros(N,1);
for k=1:N
if LES==0
Diffusion(k,1)=k^2*(1/Re)*u(k,t);
end
if LES==1
Vt_star=1+34.5*exp(-3.03*N/k);
Vt_infinite=0.31*(5-m)/(m+1)*sqrt(3-m)*Ck^(-3/2);
Ekn=u(N,t)*conj(u(N,t));
Vt=Vt_infinite*sqrt(Ekn/N)*Vt_star;
Diffusion(k,1)=k^2*(1/Re+Vt)*u(k,t);
end
for p=-N+k:N
q=k-p;
if q~=0 && p~=0
if q<0
uq=conj(u(-q,t));
else
uq=u(q,t);
end
if p<0
up=conj(u(-p,t));
else
up=u(p,t);
end
Convection(k)=Convection(k)+up*q*1i*uq;% C calculation
end
end
end
Sum(K,1)=Diffusion(K)+Convection(K);
u(K,t+1) = u(K,t)-dt*Sum(K,1);
u(1,t+1) = 1;
if max(abs(u(:,t+1)-u(:,t)))<Toler_conv
E(K) = u(K,t+1).*conj(u(K,t+1));% enegie calculation only if the error is low in off
converging_condition = 1; % Exiting condition
else
t = t+1;
Sum_ini(K) = Sum(K);
end
end
end