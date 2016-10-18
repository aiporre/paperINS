load...
('C:\Users\USER\Documents\ariel\MATLAB\Tha\plataforma2\resultados\Rtest11.mat');
M=magic(6);
Gatos=floor(size(t,2)/525)*(1:525);
Gatos(1)=1;
factor=180/pi;
M=[t(Gatos)',factor*dd4(Gatos)',factor*dd1(Gatos)',zeros(size(Gatos))'];
dlmwrite('Dphi.dat', real(M))
M=[t(Gatos)',factor*dd5(Gatos)',factor*dd2(Gatos)',zeros(size(Gatos))'];
dlmwrite('Dthe.dat',real(M))
M=[t(Gatos)',factor*dd6(Gatos)',factor*dd3(Gatos)',anglesCompleted(Gatos)'];
dlmwrite('Dpsi.dat', real(M))
