%% parte 1: Medición.
clear all
clc
h=0.01747;
t=0:h:500;
jj=59.8;
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

p_bb=p_real+0.000002533*rand(1,size(t,2)).*(sin(40*t)+cos(0.3*t))+0.000052;
q_bb=q_real+0.00000387543*rand(1,size(t,2)).*(sin(63.356*t)+cos(1.433*t))-0.0000713;
r_bb=r_real+0.0054243*rand(1,size(t,2)).*(sin(140.284*t)+cos(0.0283*t))+0.000051;
p_cc=0.0091*180*p_bb/pi+1.35;
q_cc=0.0091*180*q_bb/pi+1.35;
r_cc=0.0033*180*r_bb/pi+1.23;

axxx=gB(1,:)+0.0000002613*rand(1,size(t,2)).*(sin(10*t).*cos(0.4*t)-0.71*sin(1294.*t));
ayyy=gB(2,:)+0.0000001713*rand(1,size(t,2)).*(sin(17.8347346238430*t).*cos(28.23*t));
azzz=gB(3,:)+0.0000001813*rand(1,size(t,2)).*(sin(12.02*t).*cos(0.03*t));
ax=0.343*axxx/9.81+1.385;
ay=0.343*ayyy/9.81+1.385;
az=0.343*azzz/9.81+1.385;
%% Parte2 Variables
%Fk=[zeros(4,4),zeros(4,3);zeros(3,4),eye(3,3)];
Pk=0.1*eye(7);
%Hk=[zeros(3,4),zeros(3,3)];
%Kk=zeros(7,3);
%Qa=[0.0004,0.0004,0.0004,0.0004,0.0004,0.0004,0.0004];
%R=diag([0.987;0.99836;0.99345]);
A=zeros(4,4);
xk=[1,0,0,0,0,0,0];
p1=0;
q1=0;
r1=0;
u=0;
v=0;
w=0;
x=0;
y=0;
z=0;
%% Variables de DEBUG
aa=[];
db1=[];
db2=[];
%% Parte 3 Algoritmo iteración.
for asd=1:size(t,2)
    medicion(4)=floor(4096*p_cc(asd)/5.002);
    medicion(5)=floor(4096*q_cc(asd)/5.002);
    medicion(6)=floor(4096*r_cc(asd)/5.002);
    medicion(1)=floor(4096*ax(asd)/5.002);
    medicion(2)=floor(4096*ay(asd)/5.002);
    medicion(3)=floor(4096*az(asd)/5.002); 
    for r=1:3
      Med(r)=5.002*medicion(r)/(4096);
      Med(r)=Med(r)-1.383;
      Med(r)=9.81*Med(r)/(0.343);
    end

    for r=4:5       %%r=3;r<5;r++)
      Med(r)=5.002*medicion(r)/(4096);
      Med(r)=Med(r)-1.34;
      Med(r)=pi*Med(r)/(0.0091*180);
    end
    Med(6)=5.002*medicion(6)/(4096);
    Med(6)=Med(6)-1.226;
    Med(6)=pi*Med(6)/(0.0033*180);
    %EKF algorithm
    
    p=Med(4)-xk(5);
    q=Med(5)-xk(6);
    r=Med(6)-xk(7);
    db2=[db2 p];
    %---- s -- Hk[1][0]
    Hk(2,1)=h*sqrt(p*p+q*q+r*r)/2;%var: s
    %---I(s)--A[2)[2)
    A(3,3)=cos(Hk(2,1))+h*jj*(1-xk(1)*xk(1)-xk(2)*xk(2)-xk(3)*xk(3));
    %----Diagonal
    Fk(1,1)=A(3,3)-2*jj*h*xk(1)*xk(1);
    Fk(2,2)=A(3,3)-2*jj*h*xk(2)*xk(2);
    Fk(3,3)=A(3,3)-2*jj*h*xk(3)*xk(3);
    Fk(4,4)=A(3,3)-2*jj*h*xk(4)*xk(4);
    %----H(s)--A(1,1) y G(s)--Hk(0,0)
    if -0.0000001<Hk(2,1)<0.0000001
       A(2,2)=h/2;
       Hk(1,1)=-h/6;
    else
       A(2,2)=h*sin(Hk(2,1))/(2*Hk(2,1));
       Hk(1,1)=(h/2)*(cos(Hk(2,1))*Hk(2,1)-sin(Hk(2,1)))/(Hk(2,1)*Hk(2,1)*Hk(2,1));
    end 

    %----df/dr
    Fk(2,1)=A(2,2)*p-2*jj*h*xk(2)*xk(1);
    Fk(3,1)=A(2,2)*q-2*jj*h*xk(3)*xk(1);
    Fk(4,1)=A(2,2)*r-2*jj*h*xk(4)*xk(1);

    Fk(1,2)=-A(2,2)*p-2*jj*h*xk(1)*xk(2);
    Fk(3,2)=-A(2,2)*r-2*jj*h*xk(3)*xk(2);
    Fk(4,2)=A(2,2)*q-2*jj*h*xk(4)*xk(2);

    Fk(1,3)=-A(2,2)*q-2*jj*h*xk(1)*xk(3);
    Fk(2,3)=A(2,2)*r-2*jj*h*xk(2)*xk(3);
    Fk(4,3)=-A(2,2)*p-2*jj*h*xk(4)*xk(3);

    Fk(1,4)=-A(2,2)*r-2*jj*h*xk(1)*xk(4);
    Fk(2,4)=-A(2,2)*q-2*jj*h*xk(2)*xk(4);
    Fk(3,4)=A(2,2)*p-2*jj*h*xk(3)*xk(4);

    %%%-------df/db

    Fk(1,5)=xk(2)*p+xk(3)*q+xk(4)*r;
    Fk(1,5)=h*(A(2,2)*xk(2)+Hk(1,1)*Fk(1,5))/2;

    Fk(1,6)=Fk(1,5)*q+A(2,2)*xk(3);
    Fk(1,7)=Fk(1,5)*r+A(2,2)*xk(4);
    Fk(1,5)=Fk(1,5)*p+A(2,2)*xk(2);
    %----
    Fk(2,5)=-xk(1)*p-xk(3)*r+xk(4)*q;
    Fk(2,5)=h*(A(2,2)*xk(2)+Hk(1,1)*Fk(2,5))/2;

    Fk(2,6)=Fk(2,5)*q+A(2,2)*xk(4);
    Fk(2,7)=Fk(2,5)*r-A(2,2)*xk(3);
    Fk(2,5)=Fk(2,5)*p-A(2,2)*xk(1);
    %----
    Fk(3,5)=-xk(1)*q+xk(2)*r-xk(4)*p;
    Fk(3,5)=h*(A(2,2)*xk(3)+Hk(1,1)*Fk(3,5))/2;

    Fk(3,6)=Fk(3,5)*q-A(2,2)*xk(3);
    Fk(3,7)=Fk(3,5)*r+A(2,2)*xk(2);
    Fk(3,5)=Fk(3,5)*p-A(2,2)*xk(4);
    %----
    Fk(4,5)=-xk(1)*r-xk(2)*q+xk(3)*p;
    Fk(4,5)=h*(A(2,2)*xk(4)+Hk(1,1)*Fk(4,5))/2;

    Fk(4,6)=Fk(4,5)*q-A(2,2)*xk(2);
    Fk(4,7)=Fk(4,5)*r-A(2,2)*xk(1);
    Fk(4,5)=Fk(4,5)*p+A(2,2)*xk(3);
    %-----------------------------------------------------------------------------
    % ESTImaDO A Priori
    %-----------------------------------------------------------------------------
    A(1,1)=A(3,3)*xk(1)-A(2,2)*(p*xk(2)+q*xk(3)+r*xk(4));
    A(1,2)=A(3,3)*xk(2)-A(2,2)*(-p*xk(1)-r*xk(3)+q*xk(4));
    A(1,3)=A(3,3)*xk(3)-A(2,2)*(-q*xk(1)+r*xk(2)-p*xk(4));
    A(1,4)=A(3,3)*xk(4)-A(2,2)*(-r*xk(1)-q*xk(2)+p*xk(3));
    for i=1:4%(i=0;i<=3;i++)
      xk(i)=A(1,i);
    end
    %------------------------------------------------------------------------------
    %             HESSIANO DE y=h(xk,t)
    %------------------------------------------------------------------------------
    Hk(1,1)=-2*9.81*xk(3);
    Hk(1,2)=2*9.81*xk(4);
    Hk(1,3)=-2*9.81*xk(1);
    Hk(1,4)=2*9.81*xk(2);

    Hk(2,1)=2*9.81*xk(2);
    Hk(2,2)=2*9.81*xk(1);
    Hk(2,3)=2*9.81*xk(4);
    Hk(2,4)=2*9.81*xk(3);

    Hk(3,1)=2*9.81*xk(1);
    Hk(3,2)=-2*9.81*xk(2);
    Hk(3,3)=-2*9.81*xk(3);
    Hk(3,4)=-2*9.81*xk(4);
    %------------------------------------------------------------------------------
    %             Coovarianza A priori
    %------------------------------------------------------------------------------
%     Pk=magic(7)
%     Fk=[magic(4),[3,4,5;7,3.2,0.103;3.39,1.25,0.0023;0.1,2,65.2];zeros(3,4),eye(3)]
%     Fk*Pk*Fk'
    %Limpiamos A
    for i=1:4%(i=0;i<=3;i++)
      for j=1:4%(j=0;j<=3;i++)
         A(i,j)=0;
      end
    end
    %Actualización de P1
    for i=1:4  %(i=0;i<=3;i++)
      for n=1:4 %(p=0;p<=3;i++)
         for j=1:3 %(j=0;j<=2;i++)
           for k=1:3 %(k=0;k<=3;i++)
               A(i,n)=A(i,n)+Fk(i,k)*(Pk(k,j))*Fk(n,j);
               A(i,n)=A(i,n)+Fk(i,k)*(Pk(k,j+4))*Fk(n,j+4);
               A(i,n)=A(i,n)+Fk(i,k+4)*(Pk(k+4,j))*Fk(n,j);
               A(i,n)=A(i,n)+Fk(i,k+4)*(Pk(k+4,j+4))*Fk(n,j+4);
           end
         end
         for j=1:3 %(j=0;j<=2;i++)
            A(i,n)=A(i,n)+Fk(i,4)*Pk(4,j)*Fk(n,j);
            A(i,n)=A(i,n)+Fk(i,4)*Pk(4,j+4)*Fk(n,j+4);
         end
         for k=1:3 %(k=0;k<=2;i++)
            A(i,n)=A(i,n)+Fk(i,k)*Pk(k,4)*Fk(n,4);
            A(i,n)=A(i,n)+Fk(i,k+4)*Pk(k+4,4)*Fk(n,4);
         end
         A(i,n)=A(i,n)+Fk(i,4)*Fk(n,4)*Pk(4,4);
      end
    end
    %Copiamos A en P1;
    for i=1:4 %(i=0;i<=3;i++)
      for j=1:4 %(j=0;j<=3;i++)
         Pk(i,j)=A(i,j);
      end
    end
    %Limpiamos A
    for i=1:4 %(i=0;i<=3;i++)
      for j=1:4 %(j=0;j<=3;i++)
         A(i,j)=0;
      end
    end
    %Actualización del P2
    for i=1:4 %(i=0;i<=3;i++)
      for j=1:3 %(j=0;j<=2;i++)
         for k=1:3 %(k=0;k<=2;i++)
            A(i,j)=A(i,j)+Fk(i,k)*Pk(k,j+4)+Fk(i,k+4)*Pk(k+4,j+4);
         end
         A(i,j)=A(i,j)+Fk(i,4)*Pk(4,j+4);
      end
    end
    %Copiamos A en P2;
    for i=1:4 %(i=0;i<=3;i++)
      for j=1:3 %(j=0;j<=2;i++)
         Pk(i,j+4)=A(i,j);
      end
    end

    %Limpiamos A
    for i=1:4 %(i=0;i<=3;i++)
      for j=1:4 %(j=0;j<=3;i++)
         A(i,j)=0;
      end
    end
    %Actualización del P3
    for i=1:3 %(i=0;i<=2;i++)
      for j=1:4 %(j=0;j<=3;i++)
         for k=1:3 %(k=0;k<=2;i++)
            A(i,j)=A(i,j)+Fk(j,k)*Pk(i+4,k);
            A(i,j)=A(i,j)+Fk(j,k+4)*Pk(i+4,k+4);
         end
         A(i,j)=A(i,j)+Fk(j,4)*Pk(i+4,4);
      end
    end
    %Copiamos A en P3;
    for i=1:3 %(i=0;i<=2;i++)
      for j=1:4 %(j=0;j<=3;i++)
         Pk(i+4,j)=A(i,j);
      end
    end
    
%     Pk=(Fk*(Pk*Fk'));
    Pk(1,1)=Pk(1,1)+0.0000011;
    Pk(2,2)=Pk(2,2)+0.0000011;
    Pk(3,3)=Pk(3,3)+0.0000011;
    Pk(4,4)=Pk(4,4)+0.0000013;
    Pk(5,5)=Pk(5,5)+0.000000000000000861;
    Pk(6,6)=Pk(6,6)+0.000000000000000861;
    Pk(7,7)=Pk(7,7)+0.000000000000000861;

    %------------------------------------------------------------------------------
    %                            Ganancia de Kalman
    %------------------------------------------------------------------------------
    %Copiamos la matriz R en Fk*.
    
%    Kk=(Hk*(Pk*Hk')+R)^-1;
%    Kk=(Pk*Hk')*K;
    for i=1:3 %(i=0;i<=2;i++)
      for j=1:3 %(j=0;j<=2;j++)
         Fk(i,j)=0;
         if i==j
             Fk(i,j)=0.999892;
         end
      end
    end
    %se calcula [(dh/dr)P1(dh/dr)^T+R]_i,j      
    
    for i=1:3 %(m=0;m<=3;m++)
      for j=1:3 %(n=0;n<=2;n++)
         for k=1:4 %(s=0;s<=3;s++)
            for n=1:4 %(k=0;k<=3;k++)
                Fk(i,j)=Fk(i,j)+Hk(i,k)*Pk(k,n)*Hk(j,n);
            end
         end
      end
    end
    %Se calcula la adjunta de Fk* en Fk^(+3)
    Fk(1,4)=Fk(3,3)*Fk(2,2)-Fk(3,2)*Fk(2,3);
    Fk(2,4)=Fk(1,3)*Fk(3,2)-Fk(1,2)*Fk(3,3);
    Fk(3,4)=Fk(2,3)*Fk(1,2)-Fk(2,2)*Fk(1,3);

    Fk(1,5)=Fk(3,1)*Fk(2,3)-Fk(2,1)*Fk(3,3);
    Fk(2,5)=Fk(1,1)*Fk(3,3)-Fk(1,3)*Fk(3,1);
    Fk(3,5)=Fk(2,1)*Fk(1,3)-Fk(1,1)*Fk(2,3);

    Fk(1,6)=Fk(2,1)*Fk(3,2)-Fk(2,2)*Fk(3,1);
    Fk(2,6)=Fk(3,1)*Fk(1,2)-Fk(3,2)*Fk(1,1);
    Fk(3,6)=Fk(1,1)*Fk(2,2)-Fk(1,2)*Fk(2,1);
    
    % Calculamos el determinante
    A(1,1)=Fk(1,1)*Fk(1,4)+Fk(1,2)*Fk(1,5)+Fk(1,3)*Fk(1,6);

    % Calculamos la inversa en Fk(i,j) i,j=(1-3,1-3)
    if A(1,1)~=0        
        for i=1:3 %(i=0;i<=2;i++)
          for j=1:3 %(j=0;j<=2;j++)
             Fk(i,j)=Fk(i,j+3)/A(1,1);
          end
        end
    end 
    %Limpiamos A
    for i=1:4 %(i=0;i<=3;i++)
      for j=1:3 %(j=0;j<=3;i++)
         Kk(i,j)=0;
      end
    end
    % Calculamos finalmente Kk.
    for i=1:4 %(m=0;m<=3;m++)
      for j=1:3 %(n=0;n<=2;n++)
         for k=1:4 %(s=0;s<=3;s++)
            for n=1:3 %(k=0;k<=3;k++)
                Kk(i,j)=Kk(i,j)+Pk(i,k)*Hk(n,k)*Fk(n,j);
            end
         end
      end
    end
     
    %Limpiamos A
    for i=5:7 %(i=0;i<=3;i++)
      for j=1:3 %(j=0;j<=3;i++)
         Kk(i,j)=0;
      end
    end
    for i=1:3 %(m=0;m<=3;m++)
      for j=1:3 %(n=0;n<=2;n++)
         for k=1:4 %(s=0;s<=3;s++)
            for n=1:3 %(k=0;k<=3;k++)
                Kk(i+4,j)=Kk(i+4,j)+Pk(i+4,k)*Hk(n,k)*Fk(n,j);
            end
         end
      end
    end    

    %------------------------------------------------------------------------------
    %                            Estimado a posteori
    %------------------------------------------------------------------------------
    A(3,1)=Med(1)-2*9.81*(xk(2)*xk(4)-xk(1)*xk(3));
    A(3,2)=Med(2)-2*9.81*(xk(3)*xk(4)+xk(1)*xk(2));
    A(3,3)=Med(3)-9.81*(1-2*(xk(2)*xk(2)+xk(3)*xk(3)));
    %sqrt((2*9.81*(xk(2)*xk(4)-xk(1)*xk(3)))^2+(2*9.81*(xk(3)*xk(4)+xk(1)*xk(2)))^2+(9.81*(xk(1)*xk(1)-xk(2)*xk(2)-xk(3)*xk(3)+xk(4)*xk(4)))^2)
    db1=[db1 A(3,1)];
    for m=1:7 %(m=0;m<=6;m++)
      xk(m)=xk(m)+Kk(m,1)*A(3,1)+Kk(m,2)*A(3,2)+Kk(m,3)*A(3,3);
    end
    
    %------------------------------------------------------------------------------
    %                            Covariaza. a posteori
    %------------------------------------------------------------------------------
    %Copiamos la matriz P3 en A.
    for i=1:3 %(i=0;i<=2;i++)
      for j=1:4 %(j=0;j<=3;j++)
         A(i,j)=Pk(i+4,j);
      end
    end
    %DE AQUI EN ADELANTE NO SE CORRIGE EL ARCHIVO EN INS1.c
    %Part 1      
    for i=1:3 %(m=4;m<=6;m++)
      for j=1:4 %(n=0;n<=3;n++)
         for k=1:3 %(j=0;j<=3;j++)
            for n=1:4 %(k=0;k<=2;k++)
               A(i,j)=A(i,j)-Kk(i+4,k)*Hk(k,n)*Pk(n,j);
            end
         end
      end
    end
    %Copiamos la matriz A en P3.
    for i=1:3 %(i=0;i<=2;i++)
      for j=1:4 %(j=0;j<=3;j++)
         Pk(i+4,j)=A(i,j);
      end
    end

    %Copiamos la matriz P4 en A.
    for i=1:3 %(i=0;i<=2;i++)
      for j=1:3 %(j=0;j<=2;j++)
         A(i,j)=Pk(i+4,j+4);
      end
    end
    %Parte 2
    for i=1:3 %(m=4;m<=6;m++)
      for j=1:3 %(n=4;n<=6;n++)
         for k=1:3 %(j=0;j<=3;j++)
            for n=1:4 %(k=0;k<=2;k++)
               A(i,j)=A(i,j)-Kk(i+4,k)*Hk(k,n)*Pk(n,j+4);
               end
            end
         end
      end   
    %Copiamos la matriz A en P3.
    for i=1:3 %(i=0;i<=2;i++)
      for j=1:3 %(j=0;j<=2;j++)
         Pk(i+4,j+4)=A(i,j);
      end
    end
    %Copiamos la matriz P1 en A.
    for i=1:4 %(i=0;i<=3;i++)
      for j=1:4 %(j=0;j<=3;j++)
         A(i,j)=Pk(i,j);
      end
    end
    %Part 3      
    for i=1:4 %(m=0;m<=3;m++)
        for j=1:4 %(n=0;n<=3;n++)
         for k=1:3 %(j=0;j<=3;j++)
            for n=1:4 %(k=0;k<=2;k++)
               A(i,j)=A(i,j)-Kk(i,k)*Hk(k,n)*Pk(n,j);
            end
         end
        end
    end
    %Copiamos la matriz A en P1.
    for i=1:4 %(i=0;i<=3;i++)
      for j=1:4 %(j=0;j<=3;j++)
         Pk(i,j)=A(i,j);
      end   
    end
    %Copiamos la matriz P2 en A.
    for i=1:4 %(i=0;i<=3;i++)
      for j=1:3 %(j=0;j<=3;j++)
         A(i,j)=Pk(i,j+4);
      end
    end
    %Part 4      
    for i=1:4 %(m=0;m<=3;m++)
      for j=1:3 %(n=4;n<=6;n++)
         for k=1:3 %(j=0;j<=3;j++)
            for n=1:4 %(k=0;k<=2;k++)
               A(i,j)=A(i,j)-Kk(i,k)*Hk(k,n)*Pk(n,j+4);
            end
         end
      end
    end
    %Copiamos la matriz A en P3.
     for i=1:4 %(i=0;i<=3;i++)
      for j=1:3 %(j=0;j<=2;j++)
         Pk(i,j+4)=A(i,j);
      end
     end
    %FUNCIONES PARA LA ACTUALIZACION DE LA SALIDA
    %void InsActualizaSalidas(void)

    %Para econimizar espacio de la ram se reutilizan la variables, y la corres-
    %pondencia es la suguiente
    % /*
    % A(0:2,0:2)=Matriz de rotacion R(Theta)
    % Med(0:2)=a^B_x,a^b_y,a^B_zend(raw) más la gravedad en el marco fijo
    %         al cuerpo (g^B).
    % (p,q,r)^T=p,q,rend(corregida)
    % p1,q1,r1=p q y r en k-1
    % A(3,0:2)=dp/dt,dq/dt,dr/dtend(corregida)
    % A(0:2,3)=Aceleracion tangencial.
    % A(3,0:2)=Producto cruz wB(corregida) con rIMU
    % Fk(0,0:2)=Aceleracion centripeda.
    % Fk(1,0:2)=a^B_x,a^b_y,a^B_zend(corregida) más la gravedad en el marco fijo
    %            al cuerpo (g^B).
    % Fk(2,0:2)=a^I_x,a^I_y,a^I_zend(corregida) menos la gravedad en el marco Iner.
    % u,v,w = velocidades lineales respecto al marco inercial
    % x,y,z = posición del centro de gravedad del helicoptero respecto el marco
    %        inercial.
    % phi,theta,psi= angulos de Euler,
    %*/
    %MATRIZ DE ROTACION EN CUATERNIONES
    A(1,1)=1-(2*((xk(3)*xk(3))+(xk(4)*xk(4))));
    A(2,2)=1-(2*((xk(2)*xk(2))+(xk(4)*xk(4))));
    A(3,3)=1-(2*((xk(2)*xk(2))+(xk(3)*xk(3))));
    A(1,2)=2*((xk(2)*xk(3))-(xk(1)*xk(4)));
    A(1,3)=2*((xk(2)*xk(4))+(xk(1)*xk(3)));
    A(2,1)=2*((xk(2)*xk(3))+(xk(1)*xk(4)));
    A(2,3)=2*((xk(3)*xk(4))-(xk(1)*xk(2)));
    A(3,1)=2*((xk(2)*xk(4))-(xk(1)*xk(3)));
    A(3,2)=2*((xk(3)*xk(4))+(xk(1)*xk(2)));

    %ANGULOS DE POSE
    phi=atan2(A(3,2),A(3,3));
    theta=asin(-A(3,1));       
    psi=-atan2(A(2,1),A(1,1));
    
    aa=[aa psi];
    %MEDIDA CRUDA DE LA ACELERACION EN rIMU en m/s^2
    % for m=1:4 %(m=0;m<3;m++)
    %   Med(m)=Med(m)-1.5;
    % end
    %MEDIDA sin RUIDO ni BIAS DE LAs velocidades angulares de pose de la IMU
    p=Med(4)-xk(5);
    q=Med(5)-xk(6);
    r=Med(6)-xk(7);

    %Derivada de la velocidad angular respecto a cuerpo wB.QUEDE
    A(4,1)=(p-p1)/h;
    A(4,2)=(q-q1)/h;
    A(4,3)=(r-r1)/h;
    p1=p;
    q1=q;
    r1=r;

    %Aceleracion tangencial
         %rIMU=(1 2 3)
    A(1,4)=0*(A(4,2)*3)-(A(4,3)*1);
    A(2,4)=0*(A(4,1)*3)-(A(4,3)*0);
    A(3,4)=0*(A(4,1)*2)-(A(4,2)*0);

    %Aceleracion centripeda

    A(4,1)=0;%(q*3)-(r*2);
    A(4,2)=0;%(p*3)-(r*1);
    A(4,3)=0;%(p*2)-(q*1);

    Fk(1,1)=(q*A(4,3))-(r*A(3,2));
    Fk(1,2)=(p*A(4,3))-(r*A(3,1));
    Fk(1,3)=(p*A(4,2))-(q*A(3,1));



    %Aceleracion respecto al marco inercial estandar


    Fk(2,1)=Med(1)-Fk(1,1)-A(1,4);
    Fk(2,2)=Med(2)-Fk(1,2)-A(2,4);
    Fk(2,3)=Med(3)-Fk(1,3)-A(3,4);

    for i=1:3%(i=0;i<3;i++)

      for j=1%(j=0;j<1;j++)
         Fk(3,i)=0;
         for k=1:3%(k=0;k<3;k++)
            Fk(3,i)=Fk(3,i)+(A(i,k))*Fk(1,k);
            end
         end
      end
    Fk(3,3)=Fk(3,3)-9.79;

    %CALCULO DE LAS VELOCIDADES LINEALES.

    u=u+(h*Fk(3,1));
    v=v+(h*Fk(3,2));
    w=w+(h*Fk(3,3));

    %Calculo de la posicion.
    x=x+(h*u);
    y=y+(h*v);
    z=z+(h*w);
    xk
    %Fk
    %Hk
    %Pk
    %Kk
    
    
    end



