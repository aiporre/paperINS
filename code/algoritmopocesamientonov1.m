%%GENERACION DE LOS ANGULOS DE EULER Y LA ACELERACION Gravitacional
clc
display('Inicio de los vectores de medición')
tic
h=0.0025;

phi_real=angulos(1,:);
theta_real=angulos(2,:);
psi_real=angulos(3,:);
% % d_phi=0.23*(0.2*cos(0.2*t)-0.2*sin(0.2*t)).*exp(-0.2*t);
% % d_theta=0.5743*(0.24*cos(0.24*t)-0.2*sin(0.2*t)).*exp(-0.2*t);
% % d_psi=-0.53*(0.2*cos(0.2*t)-0.2*sin(0.2*t)).*exp(-0.2*t);
% phi_real=0.53*sin(0.2*t);
% theta_real=0.5743*sin(0.24*t);
% psi_real=-0.53*sin(1.2*t);
% 
% d_phi=0.53*0.2*cos(0.2*t);
% d_theta=0.5743*0.24*cos(0.24*t);
% d_psi=-0.53*1.2*cos(1.2*t);
% 
p_real=angvelo(1,:);
q_real=angvelo(2,:);
r_real=angvelo(3,:);
% Rtheta=[];
% Rtheta_real=[cos(theta_real).*cos(psi_real)];
% Rtheta_real=[Rtheta_real;sin(phi_real).*sin(theta_real).*cos(psi_real)-cos(phi_real).*sin(psi_real)];
% Rtheta_real=[Rtheta_real;cos(phi_real).*sin(theta_real).*cos(psi_real)+sin(phi_real).*sin(psi_real)];
% 
% Rtheta_real=[Rtheta_real;cos(theta_real).*sin(psi_real)];
% Rtheta_real=[Rtheta_real;sin(phi_real).*sin(theta_real).*sin(psi_real)+cos(phi_real).*cos(psi_real)];
% Rtheta_real=[Rtheta_real;cos(phi_real).*sin(theta_real).*sin(psi_real)-sin(phi_real).*cos(psi_real)];
% 
% Rtheta_real=[Rtheta_real;-sin(theta_real)];
% Rtheta_real=[Rtheta_real;sin(phi_real).*cos(theta_real)];
% Rtheta_real=[Rtheta_real;cos(phi_real).*cos(theta_real)];
% 
% gB=9.81*[-sin(theta_real);sin(phi_real).*cos(theta_real);cos(phi_real).*cos(theta_real)];
% 
p=Omg1(1,:)+0.012;
q=Omg1(2,:)+0.032;
r=Omg1(3,:)+0.035;

x_real=position(1,:);
y_real=position(2,:);
z_real=position(3,:);
% % 
u_real=velocity(1,:);
v_real=velocity(2,:);
w_real=velocity(3,:);

% % ax_real=-2*10*exp(-t).*cos(t);%+0.1*4.5*sin(2*t);
% % ay_real=2*10*exp(-t).*sin(t);%+0.1*1.5*sin(3.1*t);
% % az_real=2*10*exp(-t).*sin(t);%+0.1*2.5*sin(1.5*t);
% 
% x_real=140*sin(0.2*t);%exp(-t).*cos(t)+1;%-0.1*2.25*sin(2*t)/2;
% y_real=145*sin(0.31*t);%exp(-t).*sin(t)+1;%-0.1*1.5*sin(3.1*t)/(3.1^2);
% z_real=14.5*sin(0.0915*t);%exp(-t).*sin(t)+3;%-0.1*2.5*sin(1.5*t)/(1.5^2);
% 
% u_real=140*0.2*cos(0.2*t);%exp(-t).*(cos(t)-sin(t));%-0.1*2.25*cos(2*t);
% v_real=145*0.31*cos(0.31*t);%-exp(-t).*(cos(t)+sin(t));%-0.1*1.5*cos(3.1*t)/3.1;
% w_real=14.5*0.0915*cos(0.0915*t);%-exp(-t).*(cos(t)+sin(t));%-0.1*2.5*cos(1.5*t)/1.5;
% 
% ax_real=-140*0.2*0.2*sin(0.2*t);%-2*exp(-t).*cos(t);%+0.1*4.5*sin(2*t);
% ay_real=-145*0.31*0.31*sin(0.31*t);%2*exp(-t).*sin(t);%+0.1*1.5*sin(3.1*t);
% az_real=-14.5*0.0915*0.0915*sin(0.0915*t);%2*exp(-t).*sin(t);%+0.1*2.5*sin(1.5*t);
% 
% for i=1:size(t,2)
%     axx(i)=ax_real(i)*Rtheta_real(1,i)+ay_real(i)*Rtheta_real(2,i)+az_real(i)*Rtheta_real(3,i);
%     ayy(i)=ax_real(i)*Rtheta_real(4,i)+ay_real(i)*Rtheta_real(5,i)+az_real(i)*Rtheta_real(6,i);
%     azz(i)=ax_real(i)*Rtheta_real(7,i)+ay_real(i)*Rtheta_real(8,i)+az_real(i)*Rtheta_real(9,i);
% end
ax=Acc2(1,:);
ay=Acc2(2,:);
az=Acc2(3,:);

%%INICIACION  DE VARIABLES
%%
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
ll2=[];
ll3=[];
ll4=[];
ll5=[];
Matrix_ROTACION=[];
aux1=[];
aux2=[];
GPS_mea=[0;0;0];
pos_GPS=[];
%%
Q=[0.9e-7*eye(4),0.2*ones(4,3);zeros(3,4),0.8e-12*eye(3)];
%Q=0.17e-12*eye(7);
R=zeros(3,3);
R(1,1)=0.9891718;
R(2,2)=0.9891718;
R(3,3)=0.991718;
j=25.5;
Pk=0.1*eye(7,7);
xk=[1,0,0,0,0,0,0];
u=0;v=0;w=0;x=x_real(1);y=y_real(1);z=z_real(1);
%h=tout(2)-tout(1);
%t=tout';
% ax=a(:,1)'+0.00013*rand(1,size(t,2));
% ay=a(:,2)'+0.00013*rand(1,size(t,2));
% az=a(:,3)'+0.0003*rand(1,size(t,2));

hk=[0;0;0];
wB0=hk';
rIMU=0*[0.023;0.023;0.05];
Rtheta_est=eye(3);
tau_1=3.1;
tau_2=12.1;
k1=3*(tau_1+tau_2)/(tau_1*tau_2);
k2=9/(tau_1*tau_2);
bw_est=[0 0 0]';
tau3=1.861;
tau4=1.91;
tau5=27.501;
k3=3*(tau3*tau4+tau3*tau5+tau4*tau5)/(tau3*tau4*tau5);
k4=9*(tau3+tau4+tau5)/(tau3*tau4*tau5);
k5=27/(tau3*tau4*tau5);
p_est=zeros(3,1);
v_est=p_est;
a_est=p_est;

psi2=0;
display('Calculando los datos de navegación...')
for ii=1:size(p,2)
%% BLOQUE DE OBSERVADOR OPTIMO EKF    
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
   hk(3)=az(ii)+1-(9.81*(1-2*((xk(3)^2)+(xk(2)^2))));
   
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
   %% Bloque OBSERVADOR DE ORIENTACION
  %MEDIDA CRUDA DE LA ACELERACION EN rIMU en m/s^2
   aIMU=[ax(ii);ay(ii);az(ii)];
   
   %a_ms2=9.79*a_ms2/0.3;
   
   %MEDIDA sin RUIDO ni BIAS DE LAs velocidades angulares de pose de la IMU
   w_rads(1)=p(ii)-bw_est(1);
   w_rads(2)=q(ii)-bw_est(2);
   w_rads(3)=r(ii)-bw_est(3);
   %............
   %Corrección del filtro complemetario.
   %Rtheta=[Rtheta_real(1,ii),Rtheta_real(4,ii),Rtheta_real(7,ii);Rtheta_real(2,ii),Rtheta_real(5,ii),Rtheta_real(8,ii);Rtheta_real(3,ii),Rtheta_real(6,ii),Rtheta_real(9,ii)];
   Rtheta=Rtheta';
   Pi=0.5*(Rtheta*Rtheta_est'-Rtheta_est*Rtheta');
   alfa_R=k1*Rtheta_est'*[Pi(3,2);Pi(1,3);Pi(2,1)];
   alfa_w=-k2*alfa_R/k1;
   
   bw_est=bw_est+h*alfa_w;
   %calculo del estimado del la velocidad angular
   w_rads(1)=p(ii)-bw_est(1);
   w_rads(2)=q(ii)-bw_est(2);
   w_rads(3)=r(ii)-bw_est(3);
   
   ll5=[ll5 w_rads'];
   m=w_rads'+alfa_R;
   Sm=[0,-m(3),m(2);m(3),0,-m(1);-m(2),m(1),0];
   Rtheta_est=Rtheta_est*(eye(3)+h*Sm+(0.5*h^2)*Sm^2);%+(h^3/6)*Sm^3);
   %Calculo los ángulos de Euler
   theta=asin(-Rtheta_est(3,1));       
   if (theta~=(pi/2))||(theta~=(-pi/2))
          phi=atan(Rtheta_est(3,2)/Rtheta_est(3,3));
          if Rtheta_est(3,3)*cos(theta)<0
            phi=pi+phi;
          end
          psi=atan(Rtheta_est(2,1)/Rtheta_est(1,1));
          if Rtheta_est(1,1)*cos(theta)<0
            psi=pi+psi;
          end
   end       
   %Error Kll deber aproximarse a la identidad
   Kll=Rtheta_est*Rtheta;
   
   dd4=[dd4 phi];
   dd5=[dd5 theta];
   dd6=[dd6 psi];
   
 %% Bloque de combinación
   aI=Rtheta_est*aIMU;
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
  if  mod(t(ii),0.2)==0
      x=0.4*x+0.6*(pos_gps(1,floor(ii/80)+1));
      y=0.5*y+0.5*(pos_gps(2,floor(ii/80)+1));
      z=0.3*z+0.7*(pos_gps(3,floor(ii/80)+1));
  end
  pos=[x y z]';
  m=w_rads;
  Sm=[0,-m(3),m(2);m(3),0,-m(1);-m(2),m(1),0];
  alfa_p=k3*(pos-p_est);
  alfa_v=k4*(pos-p_est);
  alfa_a=-k5*(eye(3)+(1/k3)*Sm)*Rtheta'*(pos-p_est);
  p_est=p_est+(h*(v_est+alfa_p));
  v_est=v_est+(h*(aI-Rtheta*a_est+alfa_v));
  a_est=a_est+(h*alfa_a);
%   randal=[u v w]'+v_est;
%   u=0.5*randal(1);
%   v=0.5*randal(2);
%   w=0.5*randal(3);
  a_final=Rtheta_est*(aIMU-a_est)-[0;0;9.81];
  ll2=[ll2 p_est];
  ll3=[ll3 v_est];
  ll4=[ll4 a_final];
  
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
%  randal=xk(5:7)'+bw_est;
%  xk(5:7)=0.5*randal';
%  bw_est=0.5*randal;
% %   K
end
%%
toc
beep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
plot(t,ll1(1,:),'g')
plot(t,ll4(1,:),'m')
hold on
plot(t,ax_real)
hold off
title('ax en marco inercial')
pause

plot(t,ll1(2,:),'g')
plot(t,ll4(2,:),'m')
hold on
plot(t,ay_real,'b')
hold off
title('ay en marco inercial')
pause

plot(t,ll1(3,:),'g')
plot(t,ll4(3,:),'m')
hold on
plot(t,az_real,'b')
hold off
title('az en marco inercial')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(t,ll2(1,:),'m')
hold on
plot(t,aux2(1,:),'g')
plot(t,x_real,'b')
hold off
title('x en marco inercial')
pause

plot(t,ll2(2,:),'m')
hold on
plot(t,aux2(2,:),'g')
plot(t,y_real,'b')
hold off
title('y en marco inercial')
pause

plot(t,ll2(3,:),'m')
hold on
plot(t,aux2(3,:),'g')
plot(t,z_real,'b')
hold off
title('z en marco inercial')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(t,ll3(1,:),'m')
hold on
plot(t,aux1(1,:),'g')
plot(t,u_real,'b')
hold off
title('Vx en marco inercial')
pause

plot(t,ll3(2,:),'m')
hold on
plot(t,aux1(2,:),'g')
plot(t,v_real,'b')
hold off
title('Vy en marco inercial')
pause

plot(t,ll3(3,:),'m')
hold on
plot(t,aux1(3,:),'g')
plot(t,w_real,'b')
hold off
title('Vz en marco inercial')
pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
plot(t,ll5(1,:),'m')
hold off
title('velocida de alabeo')
pause


pause
plot(t,q,'r')
hold on
q2=q-uu6;
plot(t,q_real,'b')
plot(t,q2,'g')
plot(t,ll5(2,:),'m')
hold off
title('velocida de roll')
pause

pause
plot(t,r,'r')
hold on
r2=r-uu7;
plot(t,r_real,'b')
plot(t,r2,'g')
plot(t,ll5(3,:),'m')
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