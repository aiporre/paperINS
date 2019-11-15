% 
% 
% %%
% clc
% subplot(3,1,1)
% plot(t,ll2(1,:),'m')
% hold on
% plot(t,aux2(1,:),'g')
% plot(t,x_real,'b')
% plot(t0,pos_gps(1,:),'r')
% hold off
% title('x en marco inercial')
% legend('Estimación','Valor real')
% 
% subplot(3,1,2)
% plot(t,ll2(2,:),'m')
% hold on
% plot(t,aux2(2,:),'g')
% plot(t,y_real,'b')
% plot(t0,pos_gps(2,:),'r')
% hold off
% title('y en marco inercial')
% 
% subplot(3,1,3)
% plot(t,ll2(3,:),'m')
% hold on
% plot(t,aux2(3,:),'g')
% plot(t,z_real,'b')
% plot(t0,pos_gps(3,:),'r')
% hold off
% title('z en marco inercial')
% pause
% close
% %
% subplot(3,1,1)
% plot(t,ll3(1,:),'m')
% hold on
% plot(t,aux1(1,:),'g')
% plot(t,u_real,'b')
% hold off
% title('Vx en marco inercial')
% legend('Estimación','Valor real')
% 
% subplot(3,1,2)
% plot(t,ll3(2,:),'m')
% hold on
% plot(t,aux1(2,:),'g')
% plot(t,v_real,'b')
% hold off
% title('Vy en marco inercial')
% 
% subplot(3,1,3)
% plot(t,ll3(3,:),'m')
% hold on
% plot(t,aux1(3,:),'g')
% plot(t,w_real,'b')
% hold off
% title('Vz en marco inercial')
% pause
% close
% %
% subplot(3,1,1)
% plot(t,p,'r')
% hold on
% plot(t,p_real,'b')
% plot(t,ll5(1,:),'m')
% hold off
% title('velocida de alabeo')
% legend('Valor real','Estimación')
% 
% 
% subplot(3,1,2)
% plot(t,q,'r')
% hold on
% plot(t,q_real,'b')
% plot(t,ll5(2,:),'m')
% hold off
% title('velocida de roll')
% 
% 
% subplot(3,1,3)
% plot(t,r,'r')
% hold on
% plot(t,r_real,'b')
% plot(t,ll5(3,:),'m')
% hold off
% title('velocidad de giñeo')
% pause
% close
% 
% %
% subplot(3,1,1)
% plot(t,phi_real,'b')
% hold on
% plot(t,dd4,'m')
% title('Estimacion del algulode alabeo')
% legend('Valor real','Estimación')
% hold off
% 
% 
% subplot(3,1,2)
% plot(t,theta_real,'b')
% hold on
% plot(t,dd5,'m')
% title('Estimacion del algulode rodadura')
% hold off
% 
% subplot(3,1,3)
% plot(t,psi_real,'b')
% hold on
% plot(t,dd6,'m')
% title('Estimacion del algulode guiñeo')
% hold off

%%
clc
subplot(3,1,1)
plot(t,ll2(1,:),'b')
hold on
plot(t,x_real,'r')
%plot(t,pos_gps(1,:),'r')
hold off
title('x en marco inercial')
legend('Estimación','Valor real')

subplot(3,1,2)
plot(t,ll2(2,:),'b')
hold on
plot(t,y_real,'r')
%plot(t,pos_gps(2,:),'r')
hold off
title('y en marco inercial')

subplot(3,1,3)
plot(t,ll2(3,:),'b')
hold on
plot(t,z_real,'r')
%plot(t,pos_GPS(3,:),'r')
hold off
title('z en marco inercial')
pause
close
%%
subplot(3,1,1)
plot(t,ll3(1,:),'b')
hold on
plot(t,u_real,'r')
hold off
title('Vx en marco inercial')
legend('Estimación','Valor real')

subplot(3,1,2)
plot(t,ll3(2,:),'b')
hold on
plot(t,v_real,'r')
hold off
title('Vy en marco inercial')

subplot(3,1,3)
plot(t,ll3(3,:),'b')
hold on
plot(t,w_real,'r')
hold off
title('Vz en marco inercial')
pause
close
%%
subplot(3,1,1)
%plot(t,p,'r')
hold on
plot(t,p_real,'b')
hold on
plot(t,ll5(1,:),'m')
hold off
title('velocida de alabeo')
legend('Valor real','Estimación')


subplot(3,1,2)
%plot(t,q,'r')
hold on
plot(t,q_real,'b')
plot(t,ll5(2,:),'m')
hold off
title('velocida de roll')


subplot(3,1,3)
%plot(t,r,'r')
hold on
plot(t,r_real,'b')
plot(t,ll5(3,:),'m')
hold off
title('velocidad de giñeo')
pause
close

%%
subplot(3,1,1)
hold on
plot(t,phi_real,'b')
plot(t,dd1,'m')
title('Estimacion del algulode alabeo')
legend('Valor real','Estimación')
hold off


subplot(3,1,2)
plot(t,theta_real,'b')
hold on
plot(t,dd2,'m')
title('Estimacion del algulode rodadura')
hold off

subplot(3,1,3)
plot(t,psi_real,'b')

hold on
plot(t,dd3,'m')
title('Estimacion del algulode guiñeo')
hold off
pause
close


