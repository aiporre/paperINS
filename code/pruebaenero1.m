%%
asd=[];
%%

Q=eye(7,7);
Q=Q*0.012;
R=zeros(3,3);
R(1,1)=0.91;
R(2,2)=0.91;
R(3,3)=0.1;
Pk=0.1*eye(7,7);
xk=[1,0,0,0,0.10,0.10,0.10];
u=0;v=0;w=0;x=0;y=0;z=0;
h=0.2;
t=0:h:50;
p=(0.51*sin(t))+0.03*rand(1,size(t,2))+0.2;
q=(0.632*sin(t))+0.03*rand(1,size(t,2))+0.3;
r=(0.723*sin(t*3))+0.003*rand(1,size(t,2))+0.1;
ax=sin(t*10)+0.0003*rand(1,size(t,2))+0.2;
ay=sin(t*14)+0.0003*rand(1,size(t,2))+0.3;
az=sin(t*3)+0.0003*rand(1,size(t,2))+0.1;
hk=[0;0;0];
wB0=hk';
rIMU=[0.023;0.023;0.05];
for ii=1:size(p,2)
%El Jacobiano de dx/dt=f(x,n)
   Fk(1,2)=0.5*(xk(5)-p(ii));
   Fk(2,1)=-Fk(1,2);
   Fk(3,4)=-Fk(1,2);
   Fk(4,3)=Fk(1,2);
   
   Fk(1,3)=0.5*(xk(6)-q(ii));
   Fk(3,1)=-Fk(1,3);
   Fk(4,2)=-Fk(1,3);
   Fk(2,4)=Fk(1,3);
   
   Fk(1,4)=0.5*(xk(7)-r(ii));
   Fk(4,1)=-Fk(1,4);
   Fk(2,3)=-Fk(1,4);
   Fk(3,2)=Fk(1,4);
   
   Fk(1,5)=0.5*xk(2);
   Fk(1,6)=0.5*xk(3);
   Fk(1,7)=0.5*xk(4);
   
   Fk(2,5)=-0.5*xk(1);
   Fk(2,6)=0.5*xk(4);
   Fk(2,7)=-0.5*xk(3);
   
   Fk(3,5)=-0.5*xk(4);
   Fk(3,6)=-0.5*xk(1);
   Fk(3,7)=0.5*xk(2);
   
   Fk(4,5)=0.5*xk(3);
   Fk(4,6)=-0.5*xk(2);
   Fk(4,7)=-0.5*xk(1);
   
   Fk(1,1)=0;
   Fk(2,2)=0;
   Fk(3,3)=0;
   Fk(4,4)=0;
   
   Fk(5,1)=0;
   Fk(5,2)=0;
   Fk(5,3)=0;
   Fk(5,4)=0;
   Fk(5,5)=0;
   Fk(5,6)=0;
   Fk(5,7)=0;
   
   
   Fk(6,1)=0;
   Fk(6,2)=0;
   Fk(6,3)=0;
   Fk(6,4)=0;
   Fk(6,5)=0;
   Fk(6,6)=0;
   Fk(6,7)=0;
   
   Fk(7,1)=0;
   Fk(7,2)=0;
   Fk(7,3)=0;
   Fk(7,4)=0;
   Fk(7,5)=0;
   Fk(7,6)=0;
   Fk(7,7)=0;
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
   EXP=(cos(s)*eye(4))-((sin(s)/(2*s))*Phi_delta);
   xk(1:4)=(EXP*rr)';
   xk(5:7)=xk(5:7);
   AUX=[EXP,zeros(4,3);zeros(3,4),eye(3)];
% HesHk(void)
   Hk(1,1)=2*9.79*xk(3);
   Hk(1,2)=-2*9.79*xk(4);
   Hk(1,3)=2*9.79*xk(1);
   Hk(1,4)=-2*9.79*xk(2);
   Hk(1,5)=0;
   Hk(1,6)=0;
   Hk(1,7)=0;
   
   Hk(2,1)=-2*9.79*xk(2);
   Hk(2,2)=-2*9.79*xk(1);
   Hk(2,3)=-2*9.79*xk(4);
   Hk(2,4)=-2*9.79*xk(3);
   Hk(2,5)=0;
   Hk(2,6)=0;
   Hk(2,7)=0;
   
   Hk(3,1)=-2*9.79*xk(1);
   Hk(3,2)=2*9.79*xk(2);
   Hk(3,3)=2*9.79*xk(3);
   Hk(3,4)=2*9.79*xk(4);
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
   hk(1)=ax(ii)+(19.58*(xk(2)*xk(4)-xk(1)*xk(3)));
   hk(2)=ay(ii)+(19.58*(xk(3)*xk(4)+xk(1)*xk(2)));
   hk(3)=az(ii)+(9.79*(1-2*((xk(3)^2)+(xk(2)^2))));
   xk=xk+(K*hk)';
% CovAposteori(void)
Pk=(eye(7)-(K*Hk))*Pk;

% FUNCIONES PARA LA ACTUALIZACION DE LA SALIDA
%    MATRIZ DE ROTACION EN CUATERNIONES
   Rtheta(1,1)=1-(2*((xk(3)*xk(3))+(xk(4)*xk(4))));
   Rtheta(2,2)=1-(2*((xk(2)*xk(2))+(xk(4)*xk(4))));
   Rtheta(3,3)=1-(2*((xk(2)*xk(2))+(xk(3)*xk(3))));
   Rtheta(1,2)=2*((xk(2)*xk(3))-(xk(1)*xk(4)));
   Rtheta(1,3)=2*((xk(2)*xk(4))+(xk(1)*xk(3)));
   Rtheta(2,1)=2*((xk(2)*xk(3))+(xk(1)*xk(4)));
   Rtheta(2,3)=2*((xk(3)*xk(4))-(xk(1)*xk(2)));
   Rtheta(3,1)=2*((xk(2)*xk(4))-(xk(1)*xk(3)));
   Rtheta(3,2)=2*((xk(3)*xk(4))+(xk(1)*xk(2)));
   
   %ANGULOS DE POSE
   phi=atan(Rtheta(3,2)/Rtheta(3,3));
   theta=asin(-Rtheta(3,1));       
   psi=atan(Rtheta(2,1)/Rtheta(1,1));
   
   %MEDIDA CRUDA DE LA ACELERACION EN rIMU en m/s^2
   a_ms2=[ax(ii);ay(ii);az(ii)];
   %a_ms2=9.79*a_ms2/0.3;
   
   %MEDIDA sin RUIDO ni BIAS DE LAs velocidades angulares de pose de la IMU
   w_rads(1)=p(ii)-xk(4);
   w_rads(2)=q(ii)-xk(5);
   w_rads(3)=r(ii)-xk(6);
   
   %Derivada de la velocidad angular respecto a cuerpo wB.
   dotw_b=(w_rads-wB0)/0.2;
   wB0=w_rads;
   
   %Aceleracion centripeda
   
   a_centripeda=cross(w_rads,rIMU);
   a_centripeda=cross(w_rads,a_centripeda);
   
   %Aceleracion tangencial
   
   a_tangencial=cross(dotw_b,rIMU);
   
   %Aceleracion respecto al marco inercial estandar
   aIMU=a_ms2-a_centripeda'-a_tangencial';
   aI=Rtheta*aIMU;
   aI=aI-[0;0;9.79];
   
   %CALCULO DE LAS VELOCIDADES LINEALES.
   
  u=u+(0.2*aI(1));
  v=v+(0.2*aI(2));
  w=w+(0.2*aI(3));
  
  %Calculo de la posicion.
  x=x+(0.2*u);
  y=y+(0.2*v);
  z=z+(0.2*w);
  asd=[asd hk(1)];
  xk(1:4)*(xk(1:4))'
  K
end