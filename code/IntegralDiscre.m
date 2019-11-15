h=0.01;
T=10;
t=0:h:T;
f=2*sin(t);
g=3*cos(t);
f1=2*sin(t)+(0.32*rand(1,size(t,2)));
g1=3*cos(t)+(0.0231*rand(1,size(t,2)));

% f=-2*t.^3;
% f=2.^f;
% f=f.*sin(20*t)+0.01*sin(250*t);
% hd=(1+(t.^2));
% g=1./hd;
% f=f+g;
% K=400+((log(4))^2);
% Ireal=(2^(-2*T))*(20*cos(20*T)+log(4)*sin(20*T))/(K);
% Ireal=0.05-Ireal+((1-cos(250*T))/(2500))
plot(t,f)
vk1=0;
for i=1:size(t,2)
    vk1=(vk1)+(h*f(i));
end
display(vk1)
vk2=0;
for i=2:size(t,2)
    vk2=vk2+h*f(i-1);
end
display(vk2)
vk3=0;
vk3=vk3+((h/2)*f(1));
for i=2:size(t,2)
    vk3=vk3+((h/2)*(f(i)+f(i-1)));
end
display(vk3)
