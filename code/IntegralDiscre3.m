h=0.0025;
T=100;
t=0:h:T;
f=2*sin(t)+1;
g=3*cos(t);
df=2*cos(t);
dg=-3*sin(t);
d2f=-2*sin(t);
d2g=-3*cos(t);
f1=2*cos(t)+(0.032*(rand(1,size(t,2))-0.5));
g1=-3*sin(t)+(0.0231*(rand(1,size(t,2))-0.5));

% f=-2*t.^3;
% f=2.^f;
% f=f.*sin(20*t)+0.01*sin(250*t);
% hd=(1+(t.^2));
% g=1./hd;
% f=f+g;
% K=400+((log(4))^2);
% Ireal=(2^(-2*T))*(20*cos(20*T)+log(4)*sin(20*T))/(K);
% Ireal=0.05-Ireal+((1-cos(250*T))/(2500))
subplot(1,2,1)
plot(g,f,'r')
hold on
plot(g1,f1,'g')
hold off

X1=[];
Y1=[];
vk1=1;
for i=1:size(t,2)
    vk1=vk1+(h*f1(i));
    X1=[X1 vk1];
end
vk1=3;
for i=1:size(t,2)
    vk1=vk1+(h*g1(i));
    Y1=[Y1 vk1];
end


X2=[];
Y2=[];
vk1=0;
for i=2:size(t,2)
    vk1=vk1+(h*f1(i-1));
    X2=[X2 vk1];
end

vk1=3;
for i=2:size(t,2)
    vk1=vk1+(h*g1(i-1));
    Y2=[Y2 vk1];
end


X3=[];
Y3=[];
vk1=0;
vk1=vk1+((h/2)*f1(1));
for i=2:size(t,2)
    vk1=vk1+((h/2)*(f1(i)+f1(i-1)));
    X3=[X3 vk1];
end

vk1=3;
vk1=vk1+((h/2)*g1(1));
for i=2:size(t,2)
    vk1=vk1+((h/2)*(g1(i)+g1(i-1)));
    Y3=[Y3 vk1];
end

subplot(1,2,2)
plot(g,f,'r')
hold on
plot(Y1,X1,'b')
%plot(Y2,X2,'m')
%plot(Y3,X3,'g')
hold off



