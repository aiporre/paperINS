x_posteori(1)=25;
sigma_z=0.5;
sigma_x=0.001;
R=sigma_z^2;
Q=sigma_x^2;
P_posteori=R;
n=1000;
z(1)=30+sigma_z*rand(1);%simula la medición contaminada
for k=1:n
    %prediction
    x_apriori(k+1)=x_posteori(k);
    P_apriori(k+1)=P_posteori(k)+Q;
    %correction
    K(k+1)=P_apriori(k+1)/((P_apriori(k+1))+R);
    z(k+1)=30+sigma_z*rand(1);%simula la medición contaminada
    x_posteori(k+1)=x_apriori(k+1)+K(k+1)*(z(k+1)-x_apriori(k+1));
    P_posteori(k+1)=(1-K(k+1))*(P_apriori(k+1));
end
 plot(x_posteori,'r')
 grid
 hold on;
 plot(z)
 hold off;

