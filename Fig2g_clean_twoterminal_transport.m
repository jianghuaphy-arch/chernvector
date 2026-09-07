clc
clear all
%EXÄÜÁ¿£¬tauÅÝÀû¾ØÕó£¬N¿í¶È£¬M³¤¶È£¬qÎªºáÖáµãÊý£¬LÎªH00µÄ×îÖÕÎ¬¶È
dt=0.000001;
tau0=eye(2);taux=[0,1;1,0];tauy=[0,-complex(0,1);complex(0,1),0];tauz=[1,0;0,-1];
N1=40;M=10; Len=10; times=1; L=4*N1*M;

W=0;
q=200; EX=linspace(-0.33,-0.2,q);
%¶¨ÒåºÃHOO,HO1¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª
n=(N1-10)/2;
S=eye(N1); S(1:n,1:n)=zeros(n);  S(N1-n+1:N1,N1-n+1:N1)=zeros(n);
sigma1=-j*kron(eye(M),kron(S,eye(4))); sigma2=sigma1;%-j*eye(4*M*N1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª
A=1; mz=0;  m=2; tz=-0.3; tz2=-0.1;
T01=- (-1)*tauz; Tx1=(m/2*tauz)+A/2/j*taux;  Ty1=(m/2*tauz)+A/2/j*tauy;   Tz1=(tz*tauz);
T02=- (-1)*tauz; Tx2=(m/2*tauz)+A/2/j*taux;  Ty2=(m/2*tauz)+A/2/j*tauy;   Tz2=(tz2*tauz);
%¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª
tx=0.2; V=0.0;
S0=zeros(2);
T0=[T01-V*eye(2),S0;S0,T02+V*eye(2)]+tx*(kron(taux,tau0+taux));
Tz=[Tx1,S0;S0,Tx2]+0*(kron(taux,tau0));
Tx=[Ty1,S0;S0,Tz2]+0*(kron(taux,tau0));
Ty=[Tz1,S0;S0,Ty2]+0*(kron(taux,tau0));


Q1=zeros(N1);P=zeros(N1);Q1(1,N1)=1;P(N1,1)=1;
H00=kron(eye(N1),T0)+kron(diag(ones(1,N1-1),1),Tx)+kron(diag(ones(1,N1-1),-1),Tx')+kron(Q1,Tx')+kron(P,Tx);
H01=kron(eye(N1),Ty);
Q1=zeros(M);P=zeros(M);Q1(1,M)=1;P(M,1)=1;
H00=kron(eye(M),H00)+kron(diag(ones(1,M-1),1),H01)+kron(diag(ones(1,M-1),-1),H01');%+kron(Q1,H01')+kron(P,H01);
H01=kron(eye(M*N1),Tz);
%¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª
for it=1:times
    it,1
    tic
parfor i1=1:q
    i1,
    E0=EX(i1); 
    E=E0+dt*j;
    HE=eye(L)*E0;
    GrMM=inv(HE-H00-sigma1-W*(diag(rand(1,4*N1*M)-0.5)) );
    Gr1M=GrMM;
    for k=1:Len
        GrMM=inv(HE-H00-H01*GrMM*H01'-W*(diag(rand(1,4*N1*M)-0.5)) );
        Gr1M=Gr1M*H01'*GrMM;
    end
    GrMM=inv(inv(GrMM)-sigma2);
    Gr1M= Gr1M+Gr1M*sigma2*GrMM;
    Tr1(i1,it)=-real(trace((sigma1-sigma1')*Gr1M*(sigma2-sigma2')*Gr1M'));
end
toc
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N1=10;M=10; Len=40; L=4*N1*M;
%¶¨ÒåºÃHOO,HO1¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª
S=eye(N1); sigma1=-j*kron(eye(M),kron(S,eye(4))); sigma2=sigma1;%-j*eye(4*M*N1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª
T0=[T01-V*eye(2),S0;S0,T02+V*eye(2)]+tx*(kron(taux,tau0+taux));
Tx=[Tx1,S0;S0,Tx2]+0*(kron(taux,tau0));
Tz=[Ty1,S0;S0,Tz2]+0*(kron(taux,tau0));
Ty=[Tz1,S0;S0,Ty2]+0*(kron(taux,tau0));


Q1=zeros(N1);P=zeros(N1);Q1(1,N1)=1;P(N1,1)=1;
H00=kron(eye(N1),T0)+kron(diag(ones(1,N1-1),1),Tx)+kron(diag(ones(1,N1-1),-1),Tx');%+kron(Q1,Tx')+kron(P,Tx);
H01=kron(eye(N1),Ty);
Q1=zeros(M);P=zeros(M);Q1(1,M)=1;P(M,1)=1;
H00=kron(eye(M),H00)+kron(diag(ones(1,M-1),1),H01)+kron(diag(ones(1,M-1),-1),H01');%+kron(Q1,H01')+kron(P,H01);
H01=kron(eye(M*N1),Tz);
%¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª
for it=1:times
    it, 2
parfor i1=1:q
    E0=EX(i1);  
    E=E0+dt*j;
    HE=eye(L)*E0;
    GrMM=inv(HE-H00-sigma1-W*(diag(rand(1,4*N1*M)-0.5)) );
    Gr1M=GrMM;
    for k=1:Len
        GrMM=inv(HE-H00-H01*GrMM*H01'-W*(diag(rand(1,4*N1*M)-0.5)) );
        Gr1M=Gr1M*H01'*GrMM;
    end
    GrMM=inv(inv(GrMM)-sigma2);
    Gr1M= Gr1M+Gr1M*sigma2*GrMM;
    Tr2(i1,it)=-real(trace((sigma1-sigma1')*Gr1M*(sigma2-sigma2')*Gr1M'));
end
end
subplot(1,2,1)
hold on
plot(EX,Tr1,'r');
plot(EX,Tr2,'b');
plot(EX,Tr1+Tr2,'g.');


