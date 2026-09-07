clc
clear all

N=10; M=10; 
%phi1=0.0*pi;t=1; phi=pi*0.4;
phi1=0.0*pi;t=1; phi=pi*0.25;
 tau0=eye(2);taux=[0,1;1,0];tauy=[0,-complex(0,1);complex(0,1),0];tauz=[1,0;0,-1];


A=1; mz=0;  m=2; tz=-0.3; tc=0.1; 
tch=0.0;

T01=[- (-1)*tauz,tc*ones(2,1);tc*ones(1,2),-2]; 
Tz=[(m/2*tauz)+A/2/j*taux,tch*ones(2,1);tch*ones(1,2),1];  
Tx=[(m/2*tauz)+A/2/j*tauy,tch*ones(2,1);tch*ones(1,2),0.1];   
Ty=[(tz*tauz),tch*ones(2,1);tch*ones(1,2),1];


H00=kron(eye(N),T01)+(kron(diag(ones(1,N-1),1), Tx)+kron(diag(ones(1,N-1),-1), Tx'));
H01=kron(diag(exp(j*phi1*[1:N])),Ty);
H00=kron(eye(M),H00)+kron(diag(ones(1,M-1),1), H01)+kron(diag(ones(1,M-1),-1), H01');
H01=kron(kron(diag(exp(j*phi*[1:M])),eye(N)),Tz);

AT=-kron(eye(N),kron(eye(M),diag([2,1,-2])) );

q=401; K=linspace(-1,1,q);
for i=1:q
   kx=K(i)*pi;
   H=H00+H01*exp(j*kx)+H01'*exp(-j*kx);
   [wave,ef]=eig(H);
    E(i,:)=diag(ef);
    Dos(i,:)=real(diag(wave'*AT*wave));
end
subplot(2,2,4)
hold on
plot(K,E,'Color',[0.5,0.5,0.5])
axis([-1,1,-1,1])




% color_length=linspace(0,1,1000);
% Ctest=jet(1000);C11=[];C21=[];C31=[];Colour1=[];
% C=Ctest;
% VVV=(Dos/max(max(abs(Dos)))+1)*0.5;
% 
% x=2*N*M+2;
% 
% for i1=1:20:q
%     for i2=x-80:x+10
% 
%             hold on
%             if (VVV(i1,i2))>0.5
%             Colour=0.9;
%             else
%               Colour=0.2;
%             end
%             Colour1=[Colour1 Colour];
%             C1=interp1(color_length,C(:,1),Colour);
%             C2=interp1(color_length,C(:,2),Colour);
%             C3=interp1(color_length,C(:,3),Colour);
%             plot(K(i1),E(i1,i2) ,'o','MarkerSize',2,'MarkerFaceColor',[C1 C2 C3], 'MarkerEdgeColor',[C1 C2 C3])
%             %quiver(x(i1,i2)/pi,y(i1,i2)/pi,vy(i1,i2)/2,vx(i1,i2)/2,'r')
%     end
% end
% axis([-1,1,-1,1])