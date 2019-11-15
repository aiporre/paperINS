
for(r=0;r<3;r++)
  {  
  Med[r]=3.0*medicion[r]/(65535);
  Med[r]=Med[r]-1.383;
  Med[r]=9.81*Med[r]/(0.343);
  }

for(r=3;r<5;r++)
  {
  Med[r]=3.0*medicion[r]/(65535);
  Med[r]=Med[r]-1.34;
  Med[r]=pi*Med[r]/(0.0091*180);
  }

Med[5]=3.0*medicion[5]/(65535);
Med[5]=Med[5]-1.226;
Med[5]=pi*Med[5]/(0.0033*180);
%EKF algorithm

p=Med[3]-xk[4];
q=Med[4]-xk[5];
r=Med[5]-xk[6];
%---- s -- Hk[1][0]
Hk[1][0]=0.0025*sqrt(p*p+q*q+r*r)/2;%var: s
%---I(s)--A[2)[2)
A[2][2]=cos(Hk[1][0])+0.0025*120.34*(1-xk[0]*xk[0]-xk[1]*xk[1]-xk[2]*xk[2]);
%----Diagonal
Fk[0][0]=A[2][2]-2*120.34*0.0025*xk[0]*xk[0];
Fk[1][1]=A[2][2]-2*120.34*0.0025*xk[1]*xk[1];
Fk[2][2]=A[2][2]-2*120.34*0.0025*xk[2]*xk[2];
Fk[3][3]=A[2][2]-2*120.34*0.0025*xk[3]*xk[3];
%----H(s)--A(1,1) y G(s)--Hk(0,0)
if -0.0000001<Hk[1][0]<0.0000001
   {
   A[1][1]=0.0025/2;
   Hk[0][0]=-0.0025/6;
   }
else
   {
   A[1][1]=0.0025*sin(Hk[1][0])/(2*Hk[1][0]);
   Hk[0][0]=(0.0025/2)*(cos(Hk[1][0])*Hk[1][0]-sin(Hk[1][0]))/(Hk[1][0]*Hk[1][0]*Hk[1][0]);
   } 

%----df/dr
Fk[1][0]=A[1][1]*p-2*120.34*0.0025*xk[1]*xk[0];
Fk[2][0]=A[1][1]*q-2*120.34*0.0025*xk[2]*xk[0];
Fk[3][0]=A[1][1]*r-2*120.34*0.0025*xk[3]*xk[0];

Fk[0][1]=-A[1][1]*p-2*120.34*0.0025*xk[0]*xk[1];
Fk[2][1]=-A[1][1]*r-2*120.34*0.0025*xk[2]*xk[1];
Fk[3][1]=A[1][1]*q-2*120.34*0.0025*xk[3]*xk[1];

Fk[0][2]=-A[1][1]*q-2*120.34*0.0025*xk[0]*xk[2];
Fk[1][2]=A[1][1]*r-2*120.34*0.002*xk[1]*xk[2];
Fk[3][2]=-A[1][1]*p-2*120.34*0.0025*xk[3]*xk[2];

Fk[0][3]=-A[1][1]*r-2*120.34*0.0025*xk[0]*xk[3];
Fk[1][3]=-A[1][1]*q-2*120.34*0.0025*xk[1]*xk[3];
Fk[2][3]=A[1][1]*p-2*120.34*0.0025*xk[2]*xk[3];

%%%-------df/db

Fk[0][4]=xk[1]*p+xk[2]*q+xk[3]*r;
Fk[0][4]=0.0025*(A[1][1]*xk[1]+Hk[0][0]*Fk[0][4])/2;

Fk[0][5]=Fk[0][4]*q+A[1][1]*xk[2];
Fk[0][6]=Fk[0][4]*r+A[1][1]*xk[3];
Fk[0][4]=Fk[0][4]*p+A[1][1]*xk[1];
%----
Fk[1][4]=-xk[0]*p-xk[2]*r+xk[3]*q;
Fk[1][4]=0.0025*(A[1][1]*xk[1]+Hk[0][0]*Fk[1][4])/2;

Fk[1][5]=Fk[1][4]*q+A[1][1]*xk[3];
Fk[1][6]=Fk[1][4]*r-A[1][1]*xk[2];
Fk[1][4]=Fk[1][4]*p-A[1][1]*xk[0];
%----
Fk[2][4]=-xk[0]*q+xk[1]*r-xk[3]*p;
Fk[2][4]=0.0025*(A[1][1]*xk[2]+Hk[0][0]*Fk[2][4])/2;

Fk[2][5]=Fk[2][4]*q-A[1][1]*xk[2];
Fk[2][6]=Fk[2][4]*r+A[1][1]*xk[1];
Fk[2][4]=Fk[2][4]*p-A[1][1]*xk[3];
%----
Fk[3][4]=-xk[0]*r-xk[1]*q+xk[2]*p;
Fk[3][4]=0.0025*(A[1][1]*xk[3]+Hk[0][0]*Fk[3][4])/2;

Fk[3][5]=Fk[3][4]*q-A[1][1]*xk[1];
Fk[3][6]=Fk[3][4]*r-A[1][1]*xk[0];
Fk[3][4]=Fk[3][4]*p+A[1][1]*xk[2];
%-----------------------------------------------------------------------------
% ESTImaDO A Priori
%-----------------------------------------------------------------------------
A[0][0]=A[2][2]*xk[0]-A[1][1]*(p*xk[1]+q*xk[2]+r*xk[3]);
A[0][1]=A[2][2]*xk[1]-A[1][1]*(-p*xk[0]-r*xk[2]+q*xk[3]);
A[0][2]=A[2][2]*xk[2]-A[1][1]*(-q*xk[0]+r*xk[1]-p*xk[3]);
A[0][3]=A[2][2]*xk[3]-A[1][1]*(-r*xk[0]-q*xk[1]+p*xk[2]);
for (i=0;i<4;i++){%(i=0;i<=3;i++)
  xk[i]=A[0][i];
}
%------------------------------------------------------------------------------
%             HESSIANO DE y=h(xk,t)
%------------------------------------------------------------------------------
Hk[0][0]=-2*9.81*xk[2];
Hk[0][1]=2*9.81*xk[3];
Hk[0][2]=-2*9.81*xk[0];
Hk[0][3]=2*9.81*xk[1];

Hk[1][0]=2*9.81*xk[1];
Hk[1][1]=2*9.81*xk[0];
Hk[1][2]=2*9.81*xk[3];
Hk[1][3]=2*9.81*xk[2];

Hk[2][0]=2*9.81*xk[0];
Hk[2][1]=-2*9.81*xk[1];
Hk[2][2]=-2*9.81*xk[2];
Hk[2][3]=-2*9.81*xk[3];
%------------------------------------------------------------------------------
%             Coovarianza A priori
%------------------------------------------------------------------------------
%     Pk=magic[6]
%     Fk=[magic[3],[3,4,5;7,3.2,0.103;3.39,1.25,0.0023;0.1,2,65.2];zeros[2][3],eye[2]]
%     Fk*Pk*Fk'
%Limpiamos A
for (i=0;i<4;i++){%(i=0;i<=3;i++)
  for (j=0;j<4;j++){%(j=0;j<=3;i++)
     A[i][j]=0;
  }
}
%Actualización de P1
for (i=0;i<4;i++){  %(i=0;i<=3;i++)
  for (n=0;n<4;n++){ %(p=0;p<=3;i++)
     for (j=0;j<3;j++){ %(j=0;j<=2;i++)
       for (k=0;k<3;k++){ %(k=0;k<=3;i++)
           A[i][n]=A[i][n]+Fk[i][k]*(Pk[k][j])*Fk(n,j);
           A[i][n]=A[i][n]+Fk[i][k]*(Pk[k][j+4])*Fk[n][j+4];
           A[i][n]=A[i][n]+Fk[i][k+4]*(Pk[k+4][j])*Fk[n][j];
           A[i][n]=A[i][n]+Fk[i][k+4]*(Pk[k+4][j+4])*Fk[n][j+4];
       }
     }
     for (j=0;j<3;j++){ %(j=0;j<=2;i++)
        A[i][n]=A[i][n]+Fk[i][3]*Pk[3][j]*Fk[n][j];
        A[i][n]=A[i][n]+Fk[i][3]*Pk[3][j+4]*Fk[n][j+4];
     }
     for (k=0;k<3;k++){ %(k=0;k<=2;i++)
        A[i][n]=A[i][n]+Fk[i][k]*Pk[k][3]*Fk[n][3];
        A[i][n]=A[i][n]+Fk[i][k+4]*Pk[k+4][4)*Fk[n][3];
     }
     A[i][n]=A[i][n]+Fk[i][3]*Fk[n][3]*Pk[3][3];
  }
}
%Copiamos A en P1;
for (i=0;i<4;i++){ %(i=0;i<=3;i++)
  for (j=0;j<4;j++){ %(j=0;j<=3;i++)
     Pk[i][j]=A[i][j];
  }
}
%Limpiamos A
for (i=0;i<4;i++){ %(i=0;i<=3;i++)
  for (j=0;j<4;j++){ %(j=0;j<=3;i++)
     A[i][j]=0;
  }
}
%Actualización del P2
for (i=0;i<4;i++){ %(i=0;i<=3;i++)
  for (j=0;j<3;j++){ %(j=0;j<=2;i++)
     for (k=0;k<3;k++){ %(k=0;k<=2;i++)
        A[i][j]=A[i][j]+Fk[i][k]*Pk[k][j+4]+Fk[i][k+4]*Pk[k+4][j+4];
     }
     A[i][j]=A[i][j]+Fk[i][3]*Pk[3][j+4];
  }
}
%Copiamos A en P2;
for (i=0;i<4;i++){ %(i=0;i<=3;i++)
  for (j=0;j<3;j++){ %(j=0;j<=2;i++)
     Pk[i][j+4]=A[i][j];
  }
}

%Limpiamos A
for (i=0;i<4;i++){ %(i=0;i<=3;i++)
  for (j=0;j<4;j++){ %(j=0;j<=3;i++)
     A[i][j]=0;
  }
}
%Actualización del P3
for (i=0;i<3;i++){ %(i=0;i<=2;i++)
  for (j=0;j<4;j++){ %(j=0;j<=3;i++)
     for (k=0;k<3;k++){ %(k=0;k<=2;i++)
        A[i][j]=A[i][j]+Fk[j][k]*Pk[i+4][k];
        A[i][j]=A[i][j]+Fk[j][k+4]*Pk[i+4][k+4];
     }
     A[i][j]=A[i][j]+Fk[j][3]*Pk[i+4][4);
  }
}
%Copiamos A en P3;
for (i=0;i<3;i++){ %(i=0;i<=2;i++)
  for (j=0;j<4;j++){ %(j=0;j<=3;i++)
     Pk[i+4][j]=A[i][j];
  }
}

%     Pk=(Fk*(Pk*Fk'));
Pk[0][0]=Pk[0][0]+0.0000011;
Pk[1][1]=Pk[1][1]+0.0000011;
Pk[2][2]=Pk[2][2]+0.0000011;
Pk[3][3]=Pk[3][3]+0.0000013;
Pk[4][4]=Pk[4][4]+0.000000000000000861;
Pk[5][5]=Pk[5][5]+0.000000000000000861;
Pk[6][6]=Pk[6][6]+0.000000000000000861;

%------------------------------------------------------------------------------
%                            Ganancia de Kalman
%------------------------------------------------------------------------------
%Copiamos la matriz R en Fk*.

%    Kk=(Hk*(Pk*Hk')+R)^-1;
%    Kk=(Pk*Hk')*K;
for (i=0;i<3;i++){ %(i=0;i<=2;i++)
  for (j=0;j<3;j++){ %(j=0;j<=2;j++)
     Fk[i][j]=0;
     if i==j
         Fk[i][j]=0.999892;
     }
  }
}
%se calcula [(dh/dr)P1(dh/dr)^T+R]_i,j      

for (i=0;i<3;i++){ %(m=0;m<=3;m++)
  for (j=0;j<3;j++){ %(n=0;n<=2;n++)
     for (k=0;k<4;k++){ %(s=0;s<=3;s++)
        for (n=0;n<4;n++){ %(k=0;k<=3;k++)
            Fk[i][j]=Fk[i][j]+Hk[i][k]*Pk[k][n]*Hk[j][n];
        }
     }
  }
}
%Se calcula la adjunta de Fk* en Fk^(+3)
Fk[0][3]=Fk[2][2]*Fk[1][1]-Fk[2][1]*Fk[1][2];
Fk[1][3]=Fk[0][2]*Fk[2][1]-Fk[0][1]*Fk[2][2];
Fk[2][3]=Fk[1][2]*Fk[0][1]-Fk[1][1]*Fk[0][2];

Fk[0][4]=Fk[2][0]*Fk[1][2]-Fk[1][0]*Fk[2][2];
Fk[1][4]=Fk[0][0]*Fk[2][2]-Fk[0][2]*Fk[2][0];
Fk[2][4]=Fk[1][0]*Fk[0][2]-Fk[0][0]*Fk[1][2];

Fk[0][5]=Fk[1][0]*Fk[2][1]-Fk[1][1]*Fk[2][0];
Fk[1][5]=Fk[2][0]*Fk[0][1]-Fk[2][1]*Fk[0][0];
Fk[2][5]=Fk[0][0]*Fk[1][1]-Fk[0][1]*Fk[1][0];

% Calculamos el determinante
A[0][0]=Fk[0][0]*Fk[0][3]+Fk[0][1]*Fk[0][4]+Fk[0][2]*Fk[0][5];

% Calculamos la inversa en Fk[i][j] i,j=(1-3,1-3)
if A[0][0]~=0        
    for (i=0;i<3;i++){ %(i=0;i<=2;i++)
      for (j=0;j<3;j++){ %(j=0;j<=2;j++)
         Fk[i][j]=Fk[i][j+3)/A[0][0];
      }
    }
} 
%Limpiamos A
for (i=0;i<4;i++){ %(i=0;i<=3;i++)
  for (j=0;j<3;j++){ %(j=0;j<=3;i++)
     Kk[i][j]=0;
  }
}
% Calculamos finalmente Kk.
for (i=0;i<4;i++){ %(m=0;m<=3;m++)
  for (j=0;j<3;j++){ %(n=0;n<=2;n++)
     for (k=0;k<4;k++){ %(s=0;s<=3;s++)
        for (n=0;n<3;n++){ %(k=0;k<=3;k++)
            Kk[i][j]=Kk[i][j]+Pk[i][k]*Hk[n][k]*Fk[n][j];
        }
     }
  }
}

%Limpiamos A
for i=5:7 %(i=0;i<=3;i++)
  for (j=0;j<3;j++){ %(j=0;j<=3;i++)
     Kk[i][j]=0;
  }
}
for (i=0;i<3;i++){ %(m=0;m<=3;m++)
  for (j=0;j<3;j++){ %(n=0;n<=2;n++)
     for (k=0;k<4;k++){ %(s=0;s<=3;s++)
        for (n=0;n<3;n++){ %(k=0;k<=3;k++)
            Kk[i+4][j]=Kk[i+4][j]+Pk[i+4][k]*Hk[n][k]*Fk[n][j];
        }
     }
  }
}    

%------------------------------------------------------------------------------
%                            Estimado a posteori
%------------------------------------------------------------------------------
A[2][0]=Med[0]-2*9.81*(xk[1]*xk[3]-xk[0]*xk[2]);
A[2][1]=Med[1]-2*9.81*(xk[2]*xk[3]+xk[0]*xk[1]);
A[2][2]=Med[2]-9.81*(1-2*(xk[1]*xk[1]+xk[2]*xk[2]));
%sqrt((2*9.81*(xk[1]*xk[3]-xk[0]*xk[2]))^2+(2*9.81*(xk[2]*xk[3]+xk[0]*xk[1]))^2+(9.81*(xk[0]*xk[0]-xk[1]*xk[1]-xk[2]*xk[2]+xk[3]*xk[3]))^2)
db1=[db1 A[2][0]];
for m=1:7 %(m=0;m<=6;m++)
  xk(m)=xk(m)+Kk(m,1)*A[2][0]+Kk(m,2)*A[2][1]+Kk(m,3)*A[2][2];
}

%------------------------------------------------------------------------------
%                            Covariaza. a posteori
%------------------------------------------------------------------------------
%Copiamos la matriz P3 en A.
for (i=0;i<3;i++){ %(i=0;i<=2;i++)
  for (j=0;j<4;j++){ %(j=0;j<=3;j++)
     A[i][j]=Pk[i+4][j];
  }
}
%DE AQUI EN ADELANTE NO SE CORRIGE EL ARCHIVO EN INS1.c
%Part 1      
for (i=0;i<3;i++){ %(m=4;m<=6;m++)
  for (j=0;j<4;j++){ %(n=0;n<=3;n++)
     for (k=0;k<3;k++){ %(j=0;j<=3;j++)
        for (n=0;n<4;n++){ %(k=0;k<=2;k++)
           A[i][j]=A[i][j]-Kk[i+4][k]*Hk[k][n]*Pk[n][j];
        }
     }
  }
}
%Copiamos la matriz A en P3.
for (i=0;i<3;i++){ %(i=0;i<=2;i++)
  for (j=0;j<4;j++){ %(j=0;j<=3;j++)
     Pk[i+4][j]=A[i][j];
  }
}

%Copiamos la matriz P4 en A.
for (i=0;i<3;i++){ %(i=0;i<=2;i++)
  for (j=0;j<3;j++){ %(j=0;j<=2;j++)
     A[i][j]=Pk[i+4][j+4];
  }
}
%Parte 2
for (i=0;i<3;i++){ %(m=4;m<=6;m++)
  for (j=0;j<3;j++){ %(n=4;n<=6;n++)
     for (k=0;k<3;k++){ %(j=0;j<=3;j++)
        for (n=0;n<4;n++){ %(k=0;k<=2;k++)
           A[i][j]=A[i][j]-Kk[i+4][k]*Hk[k][n]*Pk[n][j+4];
           }
        }
     }
  }   
%Copiamos la matriz A en P3.
for (i=0;i<3;i++){ %(i=0;i<=2;i++)
  for (j=0;j<3;j++){ %(j=0;j<=2;j++)
     Pk[i+4][j+4]=A[i][j];
  }
}
%Copiamos la matriz P1 en A.
for (i=0;i<4;i++){ %(i=0;i<=3;i++)
  for (j=0;j<4;j++){ %(j=0;j<=3;j++)
     A[i][j]=Pk[i][j];
  }
}
%Part 3      
for (i=0;i<4;i++){ %(m=0;m<=3;m++)
    for (j=0;j<4;j++){ %(n=0;n<=3;n++)
     for (k=0;k<3;k++){ %(j=0;j<=3;j++)
        for (n=0;n<4;n++){ %(k=0;k<=2;k++)
           A[i][j]=A[i][j]-Kk[i][k]*Hk[k][n]*Pk[n][j];
        }
     }
    }
}
%Copiamos la matriz A en P1.
for (i=0;i<4;i++){ %(i=0;i<=3;i++)
  for (j=0;j<4;j++){ %(j=0;j<=3;j++)
     Pk[i][j]=A[i][j];
  }   
}
%Copiamos la matriz P2 en A.
for (i=0;i<4;i++){ %(i=0;i<=3;i++)
  for (j=0;j<3;j++){ %(j=0;j<=3;j++)
     A[i][j]=Pk[i][j+4];
  }
}
%Part 4      
for (i=0;i<4;i++){ %(m=0;m<=3;m++)
  for (j=0;j<3;j++){ %(n=4;n<=6;n++)
     for (k=0;k<3;k++){ %(j=0;j<=3;j++)
        for (n=0;n<4;n++){ %(k=0;k<=2;k++)
           A[i][j]=A[i][j]-Kk[i][k]*Hk[k][n]*Pk[n][j+4];
        }
     }
  }
}
%Copiamos la matriz A en P3.
 for (i=0;i<4;i++){ %(i=0;i<=3;i++)
  for (j=0;j<3;j++){ %(j=0;j<=2;j++)
     Pk[i][j+4]=A[i][j];
  }
 }
%FUNCIONES PARA LA ACTUALIZACION DE LA SALIDA
%void InsActualizaSalidas(void)

%Para econimizar espacio de la ram se reutilizan la variables, y la corres-
%pondencia es la suguiente
% /*
% A(0:2,0:2)=Matriz de rotacion R(Theta)
% Med(0:2)=a^B_x,a^b_y,a^B_z}(raw) más la gravedad en el marco fijo
%         al cuerpo (g^B).
% (p,q,r)^T=p,q,r}(corregida)
% p1,q1,r1=p q y r en k-1
% A[2][0:2)=dp/dt,dq/dt,dr/dt}(corregida)
% A(0:2,3)=Aceleracion tangencial.
% A[2][0:2)=Producto cruz wB(corregida) con rIMU
% Fk(0,0:2)=Aceleracion centripeda.
% Fk[0][0:2)=a^B_x,a^b_y,a^B_z}(corregida) más la gravedad en el marco fijo
%            al cuerpo (g^B).
% Fk[1][0:2)=a^I_x,a^I_y,a^I_z}(corregida) menos la gravedad en el marco Iner.
% u,v,w = velocidades lineales respecto al marco inercial
% x,y,z = posición del centro de gravedad del helicoptero respecto el marco
%        inercial.
% phi,theta,psi= angulos de Euler,
%*/
%MATRIZ DE ROTACION EN CUATERNIONES
A[0][0]=1-(2*((xk[2]*xk[2])+(xk[3]*xk[3])));
A[1][1]=1-(2*((xk[1]*xk[1])+(xk[3]*xk[3])));
A[2][2]=1-(2*((xk[1]*xk[1])+(xk[2]*xk[2])));
A[0][1]=2*((xk[1]*xk[2])-(xk[0]*xk[3]));
A[0][2]=2*((xk[1]*xk[3])+(xk[0]*xk[2]));
A[1][0]=2*((xk[1]*xk[2])+(xk[0]*xk[3]));
A[1][2]=2*((xk[2]*xk[3])-(xk[0]*xk[1]));
A[2][0]=2*((xk[1]*xk[3])-(xk[0]*xk[2]));
A[2][1]=2*((xk[2]*xk[3])+(xk[0]*xk[1]));

%ANGULOS DE POSE
phi=atan2(A[2][1],A[2][2]);
theta=asin(-A[2][0]);       
psi=atan2(A[1][0],A[0][0]);

aa=[aa psi];
%MEDIDA CRUDA DE LA ACELERACION EN rIMU en m/s^2
% for m=1:4 %(m=0;m<3;m++)
%   Med(m)=Med(m)-1.5;
% }
%MEDIDA sin RUIDO ni BIAS DE LAs velocidades angulares de pose de la IMU
p=Med[3]-xk[4];
q=Med[4]-xk[5];
r=Med[5]-xk[6];

%Derivada de la velocidad angular respecto a cuerpo wB.QUEDE
A[3][0]=(p-p1)/0.0025;
A[3][1]=(q-q1)/0.0025;
A[3][2]=(r-r1)/0.0025;
p1=p;
q1=q;
r1=r;

%Aceleracion tangencial
     %rIMU=(1 2 3)
A[0][3]=0*(A[3][1]*3)-(A[3][2]*1);
A[1][3]=0*(A[3][0]*3)-(A[3][2]*0);
A[2][3]=0*(A[3][0]*2)-(A[3][1]*0);

%Aceleracion centripeda

A[3][0]=0;%(q*3)-(r*2);
A[3][1]=0;%(p*3)-(r*1);
A[3][2]=0;%(p*2)-(q*1);

Fk[0][0]=(q*A[3][2])-(r*A[2][1]);
Fk[0][1]=(p*A[3][2])-(r*A[2][0]);
Fk[0][2]=(p*A[3][1])-(q*A[2][0]);



%Aceleracion respecto al marco inercial estandar


Fk[1][0]=Med[0]-Fk[0][0]-A[0][3];
Fk[1][1]=Med[1]-Fk[0][1]-A[1][3];
Fk[1][2]=Med[2]-Fk[0][2]-A[2][3];

for (i=0;i<3;i++){%(i=0;i<3;i++)

  for j=1%(j=0;j<1;j++)
     Fk[2][i]=0;
     for (k=0;k<3;k++){%(k=0;k<3;k++)
        Fk[2][i]=Fk[2][i]+(A[i][k])*Fk[0][k];
        }
     }
  }
Fk[2][2]=Fk[2][2]-9.79;

%CALCULO DE LAS VELOCIDADES LINEALES.

u=u+(0.0025*Fk[2][0]);
v=v+(0.0025*Fk[2][1]);
w=w+(0.0025*Fk[2][2]);

%Calculo de la posicion.
x=x+(0.0025*u);
y=y+(0.0025*v);
z=z+(0.0025*w);