clc
clear all
tau0=eye(2);taux=[0,1;1,0];tauy=[0,-complex(0,1);complex(0,1),0];tauz=[1,0;0,-1];
Periodic=1;


A=1; mz=0;  m=2; tz=-0.3;tz2=-0.1;
T01=- (-1)*tauz; Tx1=(m/2*tauz)+A/2/j*taux;  Ty1=(m/2*tauz)+A/2/j*tauy;   Tz1=(tz*tauz);
T02=- (-1)*tauz; Tx2=(m/2*tauz)+A/2/j*taux;  Ty2=(m/2*tauz)+A/2/j*tauy;   Tz2=(tz2*tauz);
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
tx=0.1; V=0.0;
S0=zeros(2);
T0=[T01-V*eye(2),S0;S0,T02+V*eye(2)]+tx*(kron(taux,tau0+taux));
Ty=[Tx1,S0;S0,Tx2]+0*(kron(taux,tau0));
Tz=[Ty1,S0;S0,Tz2]+0*(kron(taux,tau0));
Tx=[Tz1,S0;S0,Ty2]+0*(kron(taux,tau0));
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
N1=10;N2=N1;
M=101;  K=linspace(0,2,M);
x=2*N1*N2+2;

Q1=zeros(N1);P=zeros(N1);Q1(1,N1)=Periodic;P(N1,1)=Periodic;
H00=kron(eye(N1),T0)+kron(diag(ones(1,N1-1),1),Tx)+kron(diag(ones(1,N1-1),-1),Tx');     %+kron(Q1,Tx')+kron(P,Tx);
H01=kron(eye(N1),Ty);
H00=kron(eye(N2),H00)+kron(diag(ones(1,N2-1),1),H01)+kron(diag(ones(1,N2-1),-1),H01');     %+kron(Q1,H01')+kron(P,H01);
H01=kron(kron(eye(N2),eye(N1)),Tz);

AT=kron(eye(N2),kron(eye(N1),kron(tauz,tau0)) );

for i=1:M
    i
    kx=K(i)*pi;
    H=H00+(H01*exp(j*kx)+(H01*exp(j*kx))');
    [wave,ef]=eig(H);
    E(i,:)=diag(ef);
    Dos(i,:)=real(diag(wave'*AT*wave));
    % Dos(i,:)=sum(abs(wave).^4);
    T(:,i)=abs(wave(:,x)).^2;
    T1(:,i)=abs(wave(:,x+1)).^2;
end
wave=reshape(sum(reshape(T,4,N1*N1,M)),N1,N1,M);
wave1=reshape(sum(reshape(T1,4,N1*N1,M)),N1,N1,M);

plot(K,E,'r')

c=47
subplot(2,2,2)
hold on
plot(K,E,'r')%'Color',[0.5,0.5,0.5]);


color_length=linspace(0,1,1000);
Ctest=jet(1000);C11=[];C21=[];C31=[];Colour1=[];
C=Ctest;
VVV=(Dos/max(max(abs(Dos)))+1)/2;

for i1=1:20:M
    for i2=x-30:x+30
        if (VVV(i1,i2))<0.5
            hold on
            Colour=(VVV(i1,i2));
            Colour1=[Colour1 Colour];
            C1=interp1(color_length,C(:,1),Colour);
            C2=interp1(color_length,C(:,2),Colour);
            C3=interp1(color_length,C(:,3),Colour);
            plot(K(i1),E(i1,i2) ,'o','MarkerSize',2,'MarkerFaceColor',[C1 C2 C3], 'MarkerEdgeColor',[C1 C2 C3])
            %quiver(x(i1,i2)/pi,y(i1,i2)/pi,vy(i1,i2)/2,vx(i1,i2)/2,'r')
        end
    end
end
axis([0,2,-0.33,-0.2])

subplot(2,2,3)
mesh(wave1(:,:,6));
subplot(2,2,4)
mesh(wave1(:,:,50+4));
axis tight



