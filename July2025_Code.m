clc
clear all
Gr=10;Gc=10;R=4;Ec=0.5;M=2;Da=2;Q=4;Nt=0.5;Nb=0.8;phi=0.01;
Pr=6.8;psi=2;t=1;beta=0.4;h=0.1;k=0.001;Kc=2;Sc=0.22;
rhof=997.1;rhos=1.05e4;kappaf=0.613;
kappas=429;cpf=4179;cps=235;
betatf=2.1e-4;betacf=2.982e-4;betats=1.89e-5;betacs=1.07e-5;phi=0.01;
sigmaf=5.5e-6;sigmas=6.3e7;sigma=sigmaf+sigmas
b1=((1/((1-phi)^2.5))/(1-phi+phi*rhos/rhof));
%b1=3.0939;
b2=((1+(3*phi*(sigma-1)/((sigma+2)-(sigma-1)*phi)))/(1-phi+phi*(rhos/rhof)));
%b2=3.0939;
b3=((1-phi+phi*(rhos*betats/rhof*betatf))/(1-phi+phi*rhos/rhof));
%b3=2.8688;
b4=((1-phi+phi*(rhos*betacs/rhof*betacf))/(1-phi+phi*rhos/rhof));
%b4=466.6447;
b5=((rhof*cpf)/((1-phi)*rhof*cpf+phi*rhos*cps));
%b5=466.6447;
b6=((kappas+2*kappaf-2*phi*(kappaf-kappas))/((kappas+2*kappaf+phi*(kappaf-kappas))));
%b6=0.8817;
b7=(1/((1-phi)*rhof*cpf+phi*rhos*cps)); 
%b7=2.8537;
b8=((1/((1-phi)^2.5))/((1-phi)+phi*rhos*cps/(rhof*cpf)));
b9=((((kappas+2*kappaf-2*phi*(kappaf-kappas))/((kappas+2*kappaf+phi*(kappaf-kappas))))/((1-phi)*rhof*cpf+phi*rhos*cps))*((rhof*cpf)/(rhos*cps)));

Mva=[8 10 12 14];
for l=1:4
    Gc=Mva(l);
    lines={'k','r','y','m'};
    r=k/(h^2);

A1=(k*psi+0.5*b1*r);
A2=(h+b1*r);
A3=(h-k*psi-b1*r-((k*h*(b2*M+b1/Da))));
%A4=(r*h)+A1;
%A5=8-(12*r)-(4*N*k);

G1=(k*psi+(0.5*b5*r*(b6+4*R/3)/Pr)); %here B1=G1
G2=(h+(b5*r*(b6+4*R/3)/Pr));
G3=(h-k*psi-(b5*r*(b6+4*R/3)/Pr)+b7*k*h*Q);

H1=(k*psi+(0.5*r/Sc)); %D1=H1
H2=(h+(r/Sc));
H3=(h-k*psi-k*h*Kc-(r/Sc));

u=zeros(1, 101);
T=zeros(1, 101);
C=zeros(1, 101);
Q1=zeros(1 ,101);
Q2=zeros(1 ,101);
Q3=zeros(1, 101);
%Q4=zeros(1, 101);
B1=zeros(1 ,101);
B2=zeros(1 ,101);
B3=zeros(1 ,101);
%b4=zeros(1, 101);

u(1)=1;
T(1)=1;
C(1)=1;
u(101)=0.0;
T(101)=0.0;
C(101)=0.0;
for j=0:100
for i=2:100
    Q1(i)=(b3*k*h*Gr*T(i)+b4*k*h*Gc*C(i));
end
for i=2:100
B1(i)=(A3*u(i)+(A1-k*psi)*(u(i-1)+u(i+1))+Q1(i));
u(i)=((B1(i)+A1*u(i+1)+(A1-k*psi)*u(i-1))/A2);
end
for i=2:100
    Q2(i)=(b8*Ec*r*((u(i+1)-u(i))^2)+b9*r*Nb*(C(i+1)+C(i-1)-C(i)));
end
for i=2:100
B2(i)=(G3*T(i)+(G1-k*psi)*(T(i-1)+T(i+1))+Q2(i));
T(i)=((B2(i)+G1*T(i+1)+(G1-k*psi)*T(i-1))/G2);
end
for i=2:100
    Q3(i)=(0.5*r*Nt*(T(i+1)+T(i-1))-r*Nt*T(i));
end
for i=2:100
B3(i)=H3*C(i)+(H1-k*psi)*(C(i-1)+C(i+1))+Q3(i);
C(i)=((B3(i)+H1*C(i+1)+(H1-k*psi)*C(i-1)+0.5*r*Nt*(T(i+1)+T(i-1))-r*Nt*T(i))/H2);
end
if(j==100)

for i = 1:6:100
%         fprintf('%0.4f\t',u(i)); 
%         fprintf('%0.4f\n',T(i));
        %fprintf('%0.4f\n',C(i));
end
%Cf=((1/((A^(2/3))*(1-phi1)^(2.5)*(1-phi2)^(2.5)))*((u(2)-u(1))/h));
Cf=((1/((psi^(2))*(1-phi)^(2.5)))*((u(2)-u(1))/h));
Nu=-((b6+(4/3)*R)*((T(2)-T(1))/h));
Sh=-(C(2)-C(1))/h;
fprintf('%0.05f\t  %0.05f\t %0.05f\t',Cf,Nu,Sh);
fprintf('\t\n');
end
end
x=[0:0.1:10];
     %fprintf('\n==========\n');
     figure(1)
     plot(x,u(1,:),lines{l},'linewidth',2);
     xlabel('\bf Distance (y)');
     ylabel('\bf Velocity (u)');
     legend('\bf Gc = 8','\bf Gc = 10','\bf Gc = 12','\bf Gc = 14')
     str1={'\bf{ Pr = 6.8; Nt = 0.5; M = 2; Gr = 10; \psi = 2; Da = 2'
         'Kc = 2; Ec = 0.5; t = 1; R = 4;  \phi = 0.01, Sc=0.22}'};
     text(1,0.3,str1)
     
     hold on
     figure(2)
     plot(x,T(1,:),lines{l},'linewidth',2);
     xlabel('\bf Distance (y)');
     ylabel('\bf Temperature (\theta)');
     legend('\bf Gc = 8','\bf Gc = 10','\bf Gc = 12','\bf Gc = 14')
     str2={'\bf{ Pr = 6.8; Nt = 0.5; M = 2; Gr = 10; \psi = 2; Da = 2'
         'Kc = 2; Ec = 0.5; t = 1; R = 4;  \phi = 0.01, Sc=0.22}'};
     text(1,0.3,str2)
     
     hold on 
     figure(3)
     plot(x,C(1,:),lines{l},'linewidth',2)
     xlabel('\bf  Distance (y)');
     ylabel('\bf Concentration (C)');
     legend('\bf Gc = 8','\bf Gc = 10','\bf Gc = 12','\bf Gc = 14')
     str3={'\bf{ Pr = 6.8; Nt = 0.5; M = 2; Gr = 10; \psi = 2; Da = 2'
         'Kc = 2; Ec = 0.5; t = 1; R = 4;  \phi = 0.01, Sc=0.22}'};
     text(1,0.3,str3)
     
hold on
end
hold on
