clc;close all;
%function initialization
sirvI=@(ys,yi,beta,gama,N) (beta*ys*yi)/N-(gama*yi);
sirvR=@(yi,yr,gama) (gama*yi);
sirvS=@(ys,yi,beta,ro,yr,N) (-beta*ys*yi)/N-(ro*ys);
sirvV=@(ys,ro) ro*ys;
%initialize parameter
a=0;            %starting point
b=120; %end point
h=0.01;         %step size
n=(b-a)/h;      %partition number
x=a:h:b;
N=1000;       %Total population
beta=0.09;     %rate of infection
gama=0.03;    %rate of recovery
ro=0.0005;       %rate of vaccinization
iniI=0.1*N;
iniR=0;
iniV=0;
iniS=N-iniI-iniR-iniV;
%s: suspectible
ys=zeros(1,n+1);
ys(1)=iniS;
%i: infectous
yi=zeros(1,n+1);
yi(1)=iniI;
%r: recovered
yr=zeros(1,n+1);
yr(1)=iniR;
%v: vaccinised
yv=zeros(1,n+1);
yv(1)=iniV;
%calculating by RK method
for i=1:n
    %K1
    s_k1=sirvS(ys(i),yi(i),beta,ro,yr(i),N);
    i_k1=sirvI(ys(i),yi(i),beta,gama,N);
    r_k1=sirvR(yi(i),yr(i),gama);
    v_k1=sirvV(ys(i),ro);
    %K2
    s_k2=sirvS(ys(i)+(0.5*s_k1*h),yi(i)+(0.5*i_k1*h),beta,ro,yr(i)+(0.5*r_k1*h),N);
    i_k2=sirvI(ys(i)+(0.5*s_k1*h),yi(i)+(0.5*i_k1*h),beta,gama,N);
    r_k2=sirvR(yi(i)+(0.5*i_k1*h),yr(i)+(0.5*r_k1*h),gama);
    v_k2=sirvV(ys(i)+(0.5*s_k1*h),ro);
    %K3
    s_k3=sirvS(ys(i)+(0.5*s_k2*h),yi(i)+(0.5*i_k2*h),beta,ro,yr(i)+(0.5*r_k2*h),N);
    i_k3=sirvI(ys(i)+(0.5*s_k2*h),yi(i)+(0.5*i_k2*h),beta,gama,N);
    r_k3=sirvR(yi(i)+(0.5*i_k2*h),yr(i)+(0.5*r_k2*h),gama);
    v_k3=sirvV(ys(i)+(0.5*s_k2*h),ro);
    %K4
    s_k4=sirvS(ys(i)+(s_k3*h),yi(i)+(i_k3*h),beta,ro,yr(i)+(r_k3*h),N);
    i_k4=sirvI(ys(i)+(s_k3*h),yi(i)+(i_k3*h),beta,gama,N);
    r_k4=sirvR(yi(i)+(i_k3*h),yr(i)+(r_k3*h),gama);
    v_k4=sirvV(ys(i)+(s_k3*h),ro);
    %y(i+1)
    ys(i+1)=ys(i)+((h/6)*(s_k1+(2*s_k2)+(2*s_k3)+s_k4));
    yi(i+1)=yi(i)+((h/6)*(i_k1+(2*i_k2)+(2*i_k3)+i_k4));
    yr(i+1)=yr(i)+((h/6)*(r_k1+(2*r_k2)+(2*r_k3)+r_k4));
    yv(i+1)=yv(i)+((h/6)*(v_k1+(2*v_k2)+(2*v_k3)+v_k4));
end


plot(x,ys,'b',x,yi,'r',x,yr,'g',x,yv,'black',LineWidth=2);
title("SIRV Model By setting the parameter","FontSize",10);
xlabel('Time(in days)',FontSize=15);
legend('S','I','R','V',fontsize=15);