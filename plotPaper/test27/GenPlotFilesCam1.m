clear
load...
('C:\Users\USER\Documents\ariel\MATLAB\Tha\plataforma2\resultados\Rtest27.mat');
M=magic(6);
Gatos=floor(size(t1,2)/525)*(1:525);
Gatos(1)=1;
factor=1;
M=[t1(Gatos)',DD4(Gatos)',  DD1(Gatos)',zeros(size(Gatos))'];
dlmwrite('Dphi.dat', real(M))
M=[t1(Gatos)',factor*DD5(Gatos)',factor*DD2(Gatos)',anglesCompleted(Gatos)'];
dlmwrite('Dthe.dat',real(M))
M=[t1(Gatos)',factor*DD6(Gatos)',factor*DD3(Gatos)',zeros(size(Gatos))'];
dlmwrite('Dpsi.dat', real(M))
