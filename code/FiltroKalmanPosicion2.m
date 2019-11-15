h=0.0025;
T=10;
t=0:h:T;
f=exp(-t).*cos(t)+1;%Posicion real x
g=exp(-t).*sin(t)+1;%Posicion real y
df=-exp(-t).*(cos(t)+sin(t));%velocidad real x
dg=exp(-t).*(cos(t)-sin(t));%velocidad real y
d2f=2*exp(-t).*sin(t)+(0.032*(rand(1,size(t,2))-0.5));%Aceleración contaminada X
d2g=-2*exp(-t).*cos(t)+(0.0431*(rand(1,size(t,2))-0.5));%Aceleracion contaminada Y
f1=-2*cos(t);%+(0.032*sprandsym(1,size(t,2)));
g1=-3*sin(t);%+(0.0231*sprandsym(1,size(t,2)));

subplot(1,2,1)
plot(g,f,'r')
hold on
plot(g1,f1,'g')
hold off
%% Integracion por aproximacion de Euler
X1=[];
Y1=[];
vk1=1%+0.025;
vk2=2;
for i=1:size(t,2)
    vk1=vk1+(h*d2f(i));
    vk2=vk2+(h*vk1);
    X1=[X1 vk2];
end
vk1=-1;
vk2=1;
for i=1:size(t,2)
    vk1=vk1+(h*d2g(i));
    vk2=vk2+(h*vk1);
    Y1=[Y1 vk2];
end

%% Integracion por aproximacion de diferencia hacia atrás
X2=[];
Y2=[];
vk1=1;
vk2=2;
for i=2:size(t,2)
    vk1=vk1+(h*d2f(i-1));
    vk2=vk2+(h*vk1);
    X2=[X2 vk2];
end
vk1=-1;
vk2=1;
for i=2:size(t,2)
    vk1=vk1+(h*d2g(i-1));
    vk2=vk2+(h*vk1);
    Y2=[Y2 vk2];
end

%% Integracion por aproximacion de Tustin todo en uno 
X3=[];
Y3=[];
xn1=-1;
xn2=2;
xn=d2f(1)*h^2/4;
xn2=xn1;
xn1=xn;
X3=[X3 xn];
xn=2*xn1+(d2f(2)+2*d2f(1))*h^2/4;
xn2=xn1;
xn1=xn;
X3=[X3 xn];
for n=3:size(t,2)
    xn=2*xn1-xn2+(d2f(n)+2*d2f(n-1)+d2f(n-2))*h^2/4;
    xn2=xn1;
    xn1=xn;
    X3=[X3 xn];
end
xn1=-1;
xn2=2;
xn=d2g(1)*h^2/4;
xn2=xn1;
xn1=xn;
Y3=[Y3 xn];
xn=2*xn1+(d2g(2)+2*d2g(1))*h^2/4;
xn2=xn1;
xn1=xn;
Y3=[Y3 xn];

for n=3:size(t,2)
    xn=2*xn1-xn2+(d2g(n)+2*d2g(n-1)+d2g(n-2))*h^2/4;
    xn2=xn1;
    xn1=xn;
    Y3=[Y3 xn];
end

%% Integracion por aproximacion de Tustin (paso a paso).
X4=[];
Y4=[];
xk=2;
ak=0;
vk=-1;
for i=1:size(t,2)
    ak1=ak;
    ak=d2f(i);
    vk1=vk;
    vk=vk+(h/2)*(ak+ak1);
    xk=xk+(h/2)*(vk+vk1);
    X4=[X4 xk];
end
vk=0;
xk=1;
ak=0;
vk=1;
for i=1:size(t,2)
    ak1=ak;
    ak=d2g(i);
    vk1=vk;
    vk=vk+(h/2)*(ak+ak1);
    xk=xk+(h/2)*(vk+vk1);
    Y4=[Y4 xk];
end

%% Kalman lineal
x=[2;1;0;-1;1;0;0.5;0.51;0.51];
F=[eye(3),h*eye(3),((h^2)/2)*eye(3);zeros(3,3),eye(3),h*eye(3);zeros(3,3),zeros(3,3),eye(3)];
P=0.1*eye(9);
Q=0.13*eye(9);
R=0.913*eye(3);
H=[zeros(3,3),zeros(3,3),eye(3)];
v=[d2f;d2g;zeros(size(X4))];
out=[];
das=[];
for i=1:size(v,2)
    % Predicción
    x=F*x;
    P=F*P*F'+Q;
    % Corrección
    K=P*H'*inv(H*P*H'+R);
    x=x+K*(v(:,i)-H*x);
    P=(eye(9)-K*H)*P;
    out=[out v(:,i)-x(7:9,1)];
    das=[das x(1:3,1)];
end

%% Plot
subplot(1,2,2)
plot(g,f,'r')
hold on
%plot(Y1,X1,'b')
%plot(Y2,X2,'m')
plot(Y4,X4,'g')
hold off


pause
subplot(1,1,1)
plot(t,f,'r')
hold on
plot(t,X4,'b')
plot(t,das(1,:),'g')
%plot(Y2,X2,'m')
%plot(Y3,X3,'g')
hold off

pause
subplot(1,1,1)
plot(t,g,'r')
hold on
plot(t,Y4,'b')
plot(t,das(2,:),'g')
%plot(Y2,X2,'m')
%plot(Y3,X3,'g')
hold off
% 
