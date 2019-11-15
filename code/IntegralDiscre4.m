h=0.025;
T=70*pi-0.1;
t=0:h:T;
f=2*sin(t)+1;
g=3*cos(t);
df=2*cos(t);
dg=-3*sin(t);
d2f=-2*sin(t);%+(0.32*(rand(1,size(t,2))-0.5));
d2g=-3*cos(t);%+(0.231*(rand(1,size(t,2))-0.5));
f1=-2*cos(t);%+(0.032*sprandsym(1,size(t,2)));
g1=-3*sin(t);%+(0.0231*sprandsym(1,size(t,2)));

subplot(1,2,1)
plot(g,f,'r')
hold on
plot(g1,f1,'g')
hold off

X1=[];
Y1=[];
vk1=2;
vk2=1;
for i=1:size(t,2)
    vk1=vk1+(h*d2f(i));
    vk2=vk2+(h*vk1);
    X1=[X1 vk2];
end
vk1=0.037;
vk2=3;
for i=1:size(t,2)
    vk1=vk1+(h*d2g(i));
    vk2=vk2+(h*vk1);
    Y1=[Y1 vk2];
end


X2=[];
Y2=[];
vk1=2;
vk2=1;
for i=2:size(t,2)
    vk1=vk1+(h*d2f(i-1));
    vk2=vk2+(h*vk1);
    X2=[X2 vk2];
end
vk1=0;
vk2=3;
for i=2:size(t,2)
    vk1=vk1+(h*d2g(i-1));
    vk2=vk2+(h*vk1);
    Y2=[Y2 vk2];
end


X3=[];
Y3=[];
xn1=0;
xn2=0;
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

xn1=0;
xn2=0;
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

%%
X4=[];
Y4=[];
vk=0;
xk=0;
ak=0;
vk=0;
for i=1:size(t,2)
    ak1=ak;
    ak=d2f(i);
    vk1=vk;
    vk=vk+(h/2)*(ak+ak1);
    xk=xk+(h/2)*(vk+vk1);
    X4=[X4 xk];
end
vk=0;
xk=0;
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
%%

subplot(1,2,2)
plot(g,f,'r')
hold on
%plot(Y1,X1,'b')
%plot(Y2,X2,'m')
plot(Y1,X1,'g')
hold off


pause
subplot(1,1,1)
plot(f,'r')
hold on
plot(X1,'b')
%plot(Y2,X2,'m')
%plot(Y3,X3,'g')
hold off

pause
subplot(1,1,1)
plot(g,'r')
hold on
plot(Y1,'b')
%plot(Y2,X2,'m')
%plot(Y3,X3,'g')
hold off

