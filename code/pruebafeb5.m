vx=0;
h=0.1;
ax=10*(0:0.1:100);
t=(0:0.1:100);%a(:,1)';
x1=0;
x2=0;
vx=0;
for i=1:size(ax,2)
    vx=vx+h*ax(i);
    x=2*x1-x2+h^2*ax(i);
    x1=x;
    x2=x1;
end

