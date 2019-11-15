%%GENERACION DE LOS ANGULOS DE EULER Y LA ACELERACION Gravitacional
tic
h=0.0025;
t=0:h:150;
phi_real=0.53*sin(0.2*t).*exp(-0.2*t);
theta_real=0.5743*sin(0.24*t).*exp(-0.2*t);
psi_real=-0.53*sin(0.2*t).*exp(-0.2*t);
d_phi=0.23*(0.2*cos(0.2*t)-0.2*sin(0.2*t)).*exp(-0.2*t);
d_theta=0.5743*(0.24*cos(0.24*t)-0.2*sin(0.2*t)).*exp(-0.2*t);
d_psi=-0.53*(0.2*cos(0.2*t)-0.2*sin(0.2*t)).*exp(-0.2*t);
p_real=d_phi-(sin(theta_real).*d_psi);
q_real=(d_theta.*cos(phi_real))+(d_psi.*sin(phi_real).*cos(theta_real));
r_real=(-d_theta.*sin(phi_real))+(d_psi.*cos(phi_real).*cos(theta_real));
Rtheta_real=[cos(theta_real).*cos(psi_real)];
Rtheta_real=[Rtheta_real;sin(phi_real).*sin(theta_real).*cos(psi_real)-cos(phi_real).*sin(psi_real)];
Rtheta_real=[Rtheta_real;cos(phi_real).*sin(theta_real).*cos(psi_real)+sin(phi_real).*sin(psi_real)];

Rtheta_real=[Rtheta_real;cos(theta_real).*sin(psi_real)];
Rtheta_real=[Rtheta_real;sin(phi_real).*sin(theta_real).*sin(psi_real)+cos(phi_real).*cos(psi_real)];
Rtheta_real=[Rtheta_real;cos(phi_real).*sin(theta_real).*sin(psi_real)-sin(phi_real).*cos(psi_real)];

Rtheta_real=[Rtheta_real;-sin(theta_real)];
Rtheta_real=[Rtheta_real;sin(phi_real).*cos(theta_real)];
Rtheta_real=[Rtheta_real;cos(phi_real).*cos(theta_real)];

gB=9.81*[-sin(theta_real);sin(phi_real).*cos(theta_real);cos(phi_real).*cos(theta_real)];

p=p_real+0.002533*rand(1,size(t,2)).*(sin(40*t)+cos(0.3*t))+0.052;
q=q_real+0.00387543*rand(1,size(t,2)).*(sin(63.356*t)+cos(1.433*t))-0.0713;
r=r_real+0.0054243*rand(1,size(t,2)).*(sin(140.284*t)+cos(0.0283*t))+0.051;

ax=gB(1,:)+0.2613*rand(1,size(t,2)).*(sin(10*t).*cos(0.4*t)-0.71*sin(1294.*t));
ay=gB(2,:)+0.1713*rand(1,size(t,2)).*(sin(17.8347346238430*t).*cos(28.23*t));
az=gB(3,:)-12+0.1813*rand(1,size(t,2)).*(sin(12.02*t).*cos(0.03*t));

%%INICIACION  DE VARIABLES
%%
clc
das1=[];
das2=[];
asd1=[];
asd2=[];
asd3=[];
uu1=[];
uu2=[];
uu3=[];
uu4=[];
uu5=[];
uu6=[];
uu7=[];
dd1=[];
dd2=[];
dd3=[];
dd4=[];
dd5=[];
dd6=[];
ll1=[];
Matrix_ROTACION=[];
aux1=[];
aux2=[];
%%
Q=[1e-6*eye(4),zeros(4,3);zeros(3,4),0.8e-14*eye(3)];
%Q=0.17e-12*eye(7);
R=zeros(3,3);
R(1,1)=0.9891718;
R(2,2)=0.9891718;
R(3,3)=0.991718;
Pk=0.1*eye(7,7);
xk=[1,0,0,0,0,0,0];
u=0;v=0;w=0;x=0;y=0;z=0;
%h=tout(2)-tout(1);
%t=tout';
% ax=a(:,1)'+0.00013*rand(1,size(t,2));
% ay=a(:,2)'+0.00013*rand(1,size(t,2));
% az=a(:,3)'+0.0003*rand(1,size(t,2));

hk=[0;0;0];
wB0=hk';
rIMU=0*[0.023;0.023;0.05];
j=184.9851592;
Rtheta_est=eye(3);
k1=3*(0.7+11)/(0.70*11);
k2=9/(0.7*11);
bw_est=[0 0 0]';
for ii=1:size(p,2)
%El Jacobiano de dx/dt=f(x,n)

   p1=p(ii)-xk(5);%Calculamos las velocidades angulares
   q1=q(ii)-xk(6);
   r1=r(ii)-xk(7);
   %Variables...
   s=h*sqrt(p1^2+q1^2+r1^2)/2;
   LAMB=(1-xk(1)^2-xk(2)^2-xk(3)^2-xk(4)^2);
   Hs=(h/2)*(sin(s)/s);
   Gs=(h/2)*(cos(s)*s-sin(s))/(s^3);
   Is=cos(s)+h*j*LAMB;
   
   AA(1)=xk(2)*p1+xk(3)*q1+xk(4)*r1;
   AA(2)=-xk(1)*p1-xk(3)*r1+xk(4)*q1;
   AA(3)=-xk(1)*q1+xk(2)*r1-xk(4)*p1;
   AA(4)=-xk(1)*r1-xk(2)*q1+xk(3)*p1;
   
   
   Fk(1,1)=Is-2*h*j*xk(1)^2;
   Fk(2,1)=-(2*j*h*xk(1)*xk(2))+(Hs*p1);
   Fk(3,1)=-(2*j*h*xk(1)*xk(3))+(Hs*q1);
   Fk(4,1)=-(2*j*h*xk(1)*xk(4))+(Hs*r1);
   Fk(5,1)=0;
   Fk(6,1)=0;
   Fk(7,1)=0;

   Fk(1,2)=-(2*j*h*xk(2)*xk(1))-(p1*Hs);
   Fk(2,2)=Is-(2*j*h*xk(2)^2);
   Fk(3,2)=-(2*j*h*xk(2)*xk(3))-(r1*Hs);
   Fk(4,2)=-(2*j*h*xk(2)*xk(4))+(q1*Hs);
   Fk(5,2)=0;
   Fk(6,2)=0;
   Fk(7,2)=0;
   
   Fk(1,3)=-(2*j*h*xk(3)*xk(1))-(q1*Hs);
   Fk(2,3)=-(2*j*h*xk(3)*xk(2))+(r1*Hs);
   Fk(3,3)=Is-(2*j*h*xk(3)^2);
   Fk(4,3)=-(2*j*h*xk(3)*xk(4))-(p1*Hs);
   Fk(5,3)=0;
   Fk(6,3)=0;
   Fk(7,3)=0;
   
   Fk(1,4)=-(2*j*h*xk(4)*xk(1))-(r1*Hs);
   Fk(2,4)=-(2*j*h*xk(4)*xk(2))-(q1*Hs);
   Fk(3,4)=-(2*j*h*xk(4)*xk(3))+(p1*Hs);
   Fk(4,4)=Is-(2*j*h*xk(4)*xk(4));
   Fk(5,4)=0;
   Fk(6,4)=0;
   Fk(7,4)=0;
   
   Fk(1,5)=h/2*(Hs*xk(1)+Gs*AA(1))*p1+Hs*xk(2);
   Fk(2,5)=h/2*(Hs*xk(2)+Gs*AA(2))*p1-Hs*xk(1);
   Fk(3,5)=h/2*(Hs*xk(3)+Gs*AA(3))*p1-Hs*xk(4);
   Fk(4,5)=h/2*(Hs*xk(4)+Gs*AA(4))*p1+Hs*xk(3);
   Fk(5,5)=1;
   Fk(6,5)=0;
   Fk(7,5)=0;
   
   Fk(1,6)=h/2*(Hs*xk(1)+Gs*AA(1))*q1+Hs*xk(3);
   Fk(2,6)=h/2*(Hs*xk(2)+Gs*AA(2))*q1+Hs*xk(4);
   Fk(3,6)=h/2*(Hs*xk(3)+Gs*AA(3))*q1-Hs*xk(3);
   Fk(4,6)=h/2*(Hs*xk(4)+Gs*AA(4))*q1-Hs*xk(2);
   Fk(5,6)=0;
   Fk(6,6)=1;
   Fk(7,6)=0;
   
   Fk(1,7)=h/2*(Hs*xk(1)+Gs*AA(1))*r1+Hs*xk(4);
   Fk(2,7)=h/2*(Hs*xk(2)+Gs*AA(2))*r1-Hs*xk(3);
   Fk(3,7)=h/2*(Hs*xk(3)+Gs*AA(3))*r1+Hs*xk(2);
   Fk(4,7)=h/2*(Hs*xk(4)+Gs*AA(4))*r1-Hs*xk(1);
   Fk(5,7)=0;
   Fk(6,7)=0;
   Fk(7,7)=1;
   
   
%EstApriori(float p,float q,float r) 
    %----LA VERsion para el uC debe ser corregida ya que la discretizacion
    %de \dot(x)=f(x) no es correcta.
   xk(1)=Is*xk(1)-Hs*AA(1);
   xk(2)=Is*xk(2)-Hs*AA(2);
   xk(3)=Is*xk(3)-Hs*AA(3);
   xk(4)=Is*xk(4)-Hs*AA(4);
   
% HesHk(void)
   Hk(1,1)=-2*9.81*xk(3);
   Hk(1,2)=2*9.81*xk(4);
   Hk(1,3)=-2*9.81*xk(1);
   Hk(1,4)=2*9.81*xk(2);
   Hk(1,5)=0;
   Hk(1,6)=0;
   Hk(1,7)=0;
   
   Hk(2,1)=2*9.81*xk(2);
   Hk(2,2)=2*9.81*xk(1);
   Hk(2,3)=2*9.81*xk(4);
   Hk(2,4)=2*9.81*xk(3);
   Hk(2,5)=0;
   Hk(2,6)=0;
   Hk(2,7)=0;
   
   Hk(3,1)=2*9.81*xk(1);
   Hk(3,2)=-2*9.81*xk(2);
   Hk(3,3)=-2*9.81*xk(3);
   Hk(3,4)=-2*9.81*xk(4);
   Hk(3,5)=0;
   Hk(3,6)=0;
   Hk(3,7)=0;

% CovApriori(void)
   Pk=(Fk*(Pk*Fk'))+Q;
   %Pk=(AUX*(Pk*AUX'))+Q;
% KalmanGain(void)
   K=(Hk*(Pk*Hk')+R)^-1;
   K=(Pk*Hk')*K;
% EstAposteori(float ax,float ay,float az)
   hk(1)=ax(ii)-(19.62*(xk(2)*xk(4)-xk(1)*xk(3)));
   hk(2)=ay(ii)-(19.62*(xk(3)*xk(4)+xk(1)*xk(2)));
   hk(3)=az(ii)-(9.81*(1-2*((xk(3)^2)+(xk(2)^2))));
   
   zk(1)=(19.62*(xk(2)*xk(4)-xk(1)*xk(3)));
   zk(2)=(19.62*(xk(3)*xk(4)+xk(1)*xk(2)));
   zk(3)=(9.81*(1-2*((xk(3)^2)+(xk(2)^2))));
   
   das2=[das2 sqrt(zk*zk')];
   xk=xk+(K*hk)';
% CovAposteori(void)
   Pk=(eye(7)-(K*Hk))*Pk;

% FUNCIONES PARA LA ACTUALIZACION DE LA SALIDA
%    MATRIZ DE ROTACION EN CUATERNIONES
   Rtheta(1,1)=xk(1)^2+xk(2)^2-xk(3)^2-xk(4)^2;
   Rtheta(2,2)=xk(1)^2-xk(2)^2+xk(3)^2-xk(4)^2;
   Rtheta(3,3)=xk(1)^2-xk(2)^2-xk(3)^2+xk(4)^2;
   Rtheta(2,1)=2*((xk(2)*xk(3))-(xk(1)*xk(4)));
   Rtheta(3,1)=2*((xk(2)*xk(4))+(xk(1)*xk(3)));
   Rtheta(1,2)=2*((xk(2)*xk(3))+(xk(1)*xk(4)));
   Rtheta(3,2)=2*((xk(3)*xk(4))-(xk(1)*xk(2)));
   Rtheta(1,3)=2*((xk(2)*xk(4))-(xk(1)*xk(3)));
   Rtheta(2,3)=2*((xk(3)*xk(4))+(xk(1)*xk(2)));
   Matrix_ROTACION=[Matrix_ROTACION Rtheta(:)];
   %ANGULOS DE POSE
   phi=atan2(Rtheta(3,2),Rtheta(3,3));
   theta=-asin(Rtheta(3,1));       
   psi=atan(Rtheta(1,2)/Rtheta(1,1));
   dd1=[dd1 phi];
   dd2=[dd2 theta];
   dd3=[dd3 psi];
   %MEDIDA CRUDA DE LA ACELERACION EN rIMU en m/s^2
   a_ms2=[ax(ii);ay(ii);az(ii)];
   %a_ms2=9.79*a_ms2/0.3;
   
   %MEDIDA sin RUIDO ni BIAS DE LAs velocidades angulares de pose de la IMU
   w_rads(1)=p(ii)-bw_est(1);
   w_rads(2)=q(ii)-bw_est(2);
   w_rads(3)=r(ii)-bw_est(3);
   %Corrección del filtro complemetario.
   %Rtheta=[Rtheta_real(1,ii),Rtheta_real(4,ii),Rtheta_real(7,ii);Rtheta_real(2,ii),Rtheta_real(5,ii),Rtheta_real(8,ii);Rtheta_real(3,ii),Rtheta_real(6,ii),Rtheta_real(9,ii)];
   Rtheta=Rtheta';
   Pi=0.5*(Rtheta_est'*Rtheta-Rtheta'*Rtheta_est);
   alfa_R=k1*[Pi(3,2);Pi(1,3);Pi(2,1)];
   alfa_w=-k2*alfa_R/k1;
   
   bw_est=bw_est+h*alfa_w;
   w_rads(1)=p(ii)-bw_est(1);
   w_rads(2)=q(ii)-bw_est(2);
   w_rads(3)=r(ii)-bw_est(3);
   
   m=w_rads'+alfa_R;
   Sm=[0,-m(3),m(2);m(3),0,-m(1);-m(2),m(1),0];
   Rtheta_est=Rtheta_est*(eye(3)+h*Sm+(0.5*h^2)*Sm^2);%+(h^3/6)*Sm^3);
   phi=atan(Rtheta_est(3,2)/Rtheta_est(3,3));
   theta=asin(-Rtheta_est(3,1));       
   psi=atan(Rtheta_est(2,1)/Rtheta_est(1,1));
   Kll=Rtheta*Rtheta_est';
   dd4=[dd4 phi];
   dd5=[dd5 theta];
   dd6=[dd6 psi];
%    uu=[uu w_rads(1)];
   %Derivada de la velocidad angular respecto a cuerpo wB.
   dotw_b=(w_rads-wB0)/h;
   wB0=w_rads;
   
   %Aceleracion centripeda
   
   a_centripeda=cross(w_rads,rIMU);
   a_centripeda=cross(w_rads,a_centripeda);
   
   %Aceleracion tangencial
   
   a_tangencial=cross(dotw_b,rIMU);
   
   %Aceleracion respecto al marco inercial estandar
   aIMU=a_ms2-a_centripeda'-a_tangencial';
   aI=Rtheta*aIMU;
   aI=aI-[0;0;9.81];
   ll1=[ll1 aI];
   %CALCULO DE LAS VELOCIDADES LINEALES.
   
  u=u+(h*aI(1));
  v=v+(h*aI(2));
  w=w+(h*aI(3));
  
  %Calculo de la posicion.
  x=x+(h*u);
  y=y+(h*v);
  z=z+(h*w);
  aux1=[aux1 [u;v;w]];
  aux2=[aux2 [x;y;z]];
  asd1=[asd1 zk(1)];
  asd2=[asd2 zk(2)];
  asd3=[asd3 zk(3)];
  uu1=[uu1 xk(1)];
  uu2=[uu2 xk(2)];
  uu3=[uu3 xk(3)];
  uu4=[uu4 xk(4)];
  uu5=[uu5 xk(5)];
  uu6=[uu6 xk(6)];
  uu7=[uu7 xk(7)];
  das1=[das1 xk(1:4)*(xk(1:4))'];
  randal=xk(5:7)'+bw_est;
  xk(5:7)=0.5*randal';
  bw_est=0.5*randal;
%   K
end
toc
plot(t,das1)
title('magnitud del cuaternion de rotacion r')
pause
plot(t,das2)
title('magnitud del la salida del filtro (la gravedad)')
pause
plot(t,p,'r')
hold on
p2=p-uu5;
plot(t,p_real,'b')
plot(t,p2,'g')
hold off
title('velocida de alabeo')
pause


pause
plot(t,q,'r')
hold on
q2=q-uu6;
plot(t,q_real,'b')
plot(t,q2,'g')
hold off
title('velocida de roll')
pause

pause
plot(t,r,'r')
hold on
r2=r-uu7;
plot(t,r_real,'b')
plot(t,r2,'g')
hold off

title('velocidad de giñeo')
pause


plot(t,ax,'r')
hold on
ax2=gB(1,:);
ay2=gB(2,:);
az2=gB(3,:);
plot(t,asd1,'g')
plot(t,ax2,'b')
hold off
title('aceleracion del gravedad en x para el marco F^B')
pause

plot(t,ay,'r')
hold on
plot(t,asd2,'g')
plot(t,ay2,'b')
hold off
title('aceleracion del gravedad en y para el marco F^B')
pause

plot(t,az,'r')
hold on
plot(t,asd3,'g')
plot(t,az2,'b')
title('aceleracion del gravedad en z para el marco F^B')
hold off
pause

plot(t,phi_real,'b')
hold on
plot(t,dd1,'g')
plot(t,dd4,'m')
title('Estimacion del algulode alabeo')
hold off
pause

plot(t,theta_real,'b')
hold on
plot(t,dd2,'g')
plot(t,dd5,'m')
title('Estimacion del algulode rodadura')
hold off
pause

plot(t,psi_real,'b')
hold on
plot(t,dd3,'g')
plot(t,dd6,'m')
title('Estimacion del algulode guiñeo')
hold off