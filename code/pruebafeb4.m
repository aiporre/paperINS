tic
%%GENERACION DE LOS ANGULOS DE EULER Y LA ACELERACION Gravitacional
h=0.0025;
t=0:h:70;
phi_real=0.23*sin(0.2*t);
theta_real=0.5743*sin(0.24*t);
psi_real=0.23*sin(0.8*t);
d_phi=0.23*0.2*cos(0.2*t);
d_theta=0.5743*0.24*cos(0.24*t);
d_psi=0.23*0.8*cos(0.8*t);
p_real=d_phi-(sin(theta_real).*d_psi);
q_real=(d_theta.*cos(phi_real))+(d_psi.*sin(phi_real).*cos(theta_real));
r_real=(-d_theta.*sin(phi_real))+(d_psi.*cos(phi_real).*cos(theta_real));
gB=9.81*[-sin(theta_real);sin(phi_real).*cos(theta_real);cos(phi_real).*cos(theta_real)];

p=p_real+0.0000533*rand(1,size(t,2))+0.052;
q=q_real+0.000087543*rand(1,size(t,2))-0.0713;
r=r_real+0.00004243*rand(1,size(t,2))+0.051;

ax=gB(1,:)+0.2*sin(t)+0.00613*rand(1,size(t,2));
ay=gB(2,:)+0.03*sin(t)+0.00713*rand(1,size(t,2));
az=gB(3,:)-10.02+0.00813*rand(1,size(t,2));

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
ll1=[];
Matrix_ROTACION=[];
%%
Q=[5e-6*eye(4),zeros(4,3);zeros(3,4),0.8e-12*eye(3)];
%Q=0.17e-12*eye(7);
R=zeros(3,3);
R(1,1)=0.891718;
R(2,2)=0.891718;
R(3,3)=0.91718;
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
j=100.9851592;
for ii=1:size(p,2)
%El Jacobiano de dx/dt=f(x,n)

   p1=p(ii)-xk(5);%Calculamos las velocidades angulares
   q1=q(ii)-xk(6);
   r1=r(ii)-xk(7);
   
   s=h*sqrt(p1^2+q1^2+r1^2)/2;
   
   gamma=cos(s)+j*h*(1-xk(1)^2-xk(2)^2-xk(3)^2-xk(4)^2);
   
   gamma2=-((h^3)*(s*cos(s)-sin(s)))/(8*s^3);
   phi_p=p1*gamma2;
   phi_q=q1*gamma2;
   phi_r=r1*gamma2;
   
   gamma3=((h^2)/4)*sin(s)/s;
   
   gamma4=h*sin(s)/(2*s);
   
   AA(1)=xk(2)*p1+xk(3)*q1+xk(4)*r1;
   AA(2)=-xk(1)*p1-xk(3)*r1+xk(4)*q1;
   AA(3)=-xk(1)*q1+xk(2)*r1-xk(4)*p1;
   AA(4)=-xk(1)*r1-xk(2)*q1+xk(3)*p1;
   
   
   Fk(1,1)=gamma-(2*j*h*xk(1)*xk(1));
   Fk(2,1)=-(2*j*h*xk(1)*xk(2))+(p1*h*sin(s)/(2*s));
   Fk(3,1)=-(2*j*h*xk(1)*xk(3))+(q1*h*sin(s)/(2*s));
   Fk(4,1)=-(2*j*h*xk(1)*xk(4))+(r1*h*sin(s)/(2*s));
   Fk(5,1)=0;
   Fk(6,1)=0;
   Fk(7,1)=0;

   Fk(1,2)=-(2*j*h*xk(2)*xk(1))-(p1*h*sin(s)/(2*s));
   Fk(2,2)=gamma-(2*j*h*xk(2)*xk(2));
   Fk(3,2)=-(2*j*h*xk(2)*xk(3))-(r1*h*sin(s)/(2*s));
   Fk(4,2)=-(2*j*h*xk(2)*xk(4))+(q1*h*sin(s)/(2*s));
   Fk(5,2)=0;
   Fk(6,2)=0;
   Fk(7,2)=0;
   
   Fk(1,3)=-(2*j*h*xk(3)*xk(1))-(q1*h*sin(s)/(2*s));
   Fk(2,3)=-(2*j*h*xk(3)*xk(2))+(r1*h*sin(s)/(2*s));
   Fk(3,3)=gamma-(2*j*h*xk(3)*xk(3));
   Fk(4,3)=-(2*j*h*xk(3)*xk(4))-(p1*h*sin(s)/(2*s));
   Fk(5,3)=0;
   Fk(6,3)=0;
   Fk(7,3)=0;
   
   Fk(1,4)=-(2*j*h*xk(4)*xk(1))-(r1*h*sin(s)/(2*s));
   Fk(2,4)=-(2*j*h*xk(4)*xk(2))-(q1*h*sin(s)/(2*s));
   Fk(3,4)=-(2*j*h*xk(4)*xk(3))+(p1*h*sin(s)/(2*s));
   Fk(4,4)=gamma-(2*j*h*xk(4)*xk(4));
   Fk(5,4)=0;
   Fk(6,4)=0;
   Fk(7,4)=0;
   
   Fk(1,5)=gamma3*xk(1)*p1+phi_p*AA(1)-gamma4*xk(2);
   Fk(2,5)=gamma3*xk(2)*p1+phi_p*AA(2)+gamma4*xk(1);
   Fk(3,5)=gamma3*xk(3)*p1+phi_p*AA(3)+gamma4*xk(4);
   Fk(4,5)=gamma3*xk(4)*p1+phi_p*AA(4)-gamma4*xk(3);
   Fk(5,5)=-1;
   Fk(6,5)=0;
   Fk(7,5)=0;
   
   
   Fk(1,6)=gamma3*xk(1)*q1+phi_q*AA(1)-gamma4*xk(3);
   Fk(2,6)=gamma3*xk(2)*q1+phi_q*AA(2)-gamma4*xk(4);
   Fk(3,6)=gamma3*xk(3)*q1+phi_q*AA(3)+gamma4*xk(1);
   Fk(4,6)=gamma3*xk(4)*q1+phi_q*AA(4)+gamma4*xk(2);
   Fk(5,6)=0;
   Fk(6,6)=-1;
   Fk(7,6)=0;
   
   Fk(1,7)=gamma3*xk(1)*r1+phi_r*AA(1)-gamma4*xk(4);
   Fk(2,7)=gamma3*xk(2)*r1+phi_r*AA(2)+gamma4*xk(3);
   Fk(3,7)=gamma3*xk(3)*r1+phi_r*AA(3)-gamma4*xk(2);
   Fk(4,7)=gamma3*xk(4)*r1+phi_r*AA(4)+gamma4*xk(1);
   Fk(5,7)=0;
   Fk(6,7)=0;
   Fk(7,7)=-1;
   
   
%EstApriori(float p,float q,float r) 
    %----LA VERsion para el uC debe ser corregida ya que la discretizacion
    %de \dot(x)=f(x) no es correcta.
   p1=p(ii)-xk(5);%Calculamos las velocidades angulares
   q1=q(ii)-xk(6);
   r1=r(ii)-xk(7);
   dphi=h*p1;
   dtheta=h*q1;
   dpsi=h*r1;
   Phi_delta=[0 dphi dtheta dpsi;-dphi  0 -dpsi dtheta; -dtheta dpsi 0 -dphi;-dpsi -dtheta dphi 0];
   rr=xk(1:4)';
   s=dot([dphi dtheta dpsi],[dphi dtheta dpsi]);
   s=0.5*sqrt(s);
   EXP=(gamma*eye(4))-((sin(s)/(2*s))*Phi_delta);
   xk(1:4)=(EXP*rr)';
%    asdfgh=xk(1:4)*xk(1:4)';
%    xk(1:4)=xk(1:4)/asdfgh;
   xk(5:7)=xk(5:7);
   AUX=[EXP,zeros(4,3);zeros(3,4),eye(3)];
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
   theta=asin(-Rtheta(3,1));       
   psi=atan2(Rtheta(1,2),-Rtheta(1,1));
   dd1=[dd1 phi];
   dd2=[dd2 theta];
   dd3=[dd3 psi];
   %MEDIDA CRUDA DE LA ACELERACION EN rIMU en m/s^2
   a_ms2=[ax(ii);ay(ii);az(ii)];
   %a_ms2=9.79*a_ms2/0.3;
   
   %MEDIDA sin RUIDO ni BIAS DE LAs velocidades angulares de pose de la IMU
   w_rads(1)=p(ii)-xk(4);
   w_rads(2)=q(ii)-xk(5);
   w_rads(3)=r(ii)-xk(6);
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
   aI=Rtheta'*aIMU;
   aI=aI+[0;0;9.81];
   ll1=[ll1 aI];
   %CALCULO DE LAS VELOCIDADES LINEALES.
   
  u=u+(h*aI(1));
  v=v+(h*aI(2));
  w=w+(h*aI(3));
  
  %Calculo de la posicion.
  x=x+(h*u);
  y=y+(h*v);
  z=z+(h*w);
  
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
title('Estimacion del algulode alabeo')
hold off
pause

plot(t,theta_real,'b')
hold on
plot(t,dd2,'g')
title('Estimacion del algulode rodadura')
hold off
pause

plot(t,psi_real,'b')
hold on
plot(t,dd3,'g')
title('Estimacion del algulode guiñeo')
hold off
pause
