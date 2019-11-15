T=10000;
h=0.001:0.001:0.5;
i1=[];
i2=[];
i3=[];
K=1+((log(4)*log(4))/400);
Ireal=(0.05-((2^(-2*T))*((cos(20*T)/20)+((log(4)/400)*sin(20*T)))))/(K);
Ireal=Ireal+((1-cos(250*T))/(2500));   

for j=1:size(h,2)
    t=0:h(j):T;
    f=-2*t;
    f=2.^f;
    f=f.*sin(20*t)+0.1*sin(250*t);
    hd=(1+(t.^2));
    g=1./hd;
    f=f+g;
    %plot(t,f)
    vk1=0;
    for i=1:size(t,2)
        vk1=(vk1)+(h(j)*f(i));
    end
    % display(vk1)
    vk2=0;
    for i=2:size(t,2)
        vk2=vk2+h(j)*f(i-1);
    end
    % display(vk2)
    vk3=0;
    vk3=vk3+((h(j)/2)*f(1));
    for i=2:size(t,2)
        vk3=vk3+((h(j)/2)*(f(i)+f(i-1)));
    end
    % display(vk3)
    i1=[i1 vk1];
    i2=[i2 vk2];
    i3=[i3 vk3];
end
a=ones(size(h,2));
a=a*MATLABIntegral;
% a=a*Ireal;
plot(h,i1,'r')
hold on
plot(h,i2,'g')
plot(h,i3,'b')
plot(h,a,'m')
hold off