h=0.00025;
T=100;
t=0:h:T;
f=-sin(t);%Posicion real x
g=-2*cos(2*t);%Posicion real y
df=-cos(t);%velocidad real x
dg=4*sin(2*t);%velocidad real y
d2f=sin(t);%+(0.032*(rand(1,size(t,2))-0.5));%Aceleración contaminada X
d2g=8*cos(2*t);%+(0.0431*(rand(1,size(t,2))-0.5));%Aceleracion contaminada Y

subplot(1,2,1)
plot(g,f,'r')
hold on
plot(g1,f1,'g')
hold off
%% Integracion por aproximacion de Tustin (paso a paso).
X4=[];
Y4=[];
xk=0;
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
xk=-8;
ak=0;
vk=0;
for i=1:size(t,2)
    ak1=ak;
    ak=d2g(i);
    vk1=vk;
    vk=vk+(h/2)*(ak+ak1);
    xk=xk+(h/2)*(vk+vk1);
    Y4=[Y4 xk];
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
%plot(t,das(1,:),'g')
%plot(Y2,X2,'m')
%plot(Y3,X3,'g')
hold off

pause
subplot(1,1,1)
plot(t,g,'r')
hold on
plot(t,Y4,'b')
%plot(t,das(2,:),'g')
%plot(Y2,X2,'m')
%plot(Y3,X3,'g')
hold off
% 
