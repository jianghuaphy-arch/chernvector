clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%中心区大小宽度N，长度L，，time无序循环次数，w无序强度
%所有端口宽度均相同 % N是电流电极宽度，NV电压电极宽度， NH高度
N=10; NV=10; NH=10; K=40; dt=0.0002; W=0; time=1;    sigma=1;

q=40;   Energy=linspace(-0.6,-0.1,q);
tau0=eye(2);taux=[0,1;1,0];tauy=[0,-complex(0,1);complex(0,1),0];tauz=[1,0;0,-1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%定义如下h01,hx1,hy1即可运行
A=1; mz=0;  m=2; tz=-0.3; tz2=-0.1;
T01=- (-1)*tauz; Tx1=(m/2*tauz)+A/2/j*taux;  Ty1=(m/2*tauz)+A/2/j*tauy;   Tz1=(tz*tauz);
T02=- (-1)*tauz; Tx2=(m/2*tauz)+A/2/j*taux;  Ty2=(m/2*tauz)+A/2/j*tauy;   Tz2=(tz2*tauz);
%————————————————————————————————————————
tx=0.1; V=0.0;
S0=zeros(2);
T0=[T01,S0;S0,T02]+tx*(kron(taux,tau0+taux));
Ty=[Tx1,S0;S0,Tx2]+0*(kron(taux,tau0));
Tz=[Ty1,S0;S0,Tz2]+0*(kron(taux,tau0));
Tx=[Tz1,S0;S0,Ty2]+0*(kron(taux,tau0));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h01=kron(eye(NH),T0)+kron(diag(ones(1,NH-1),1),Tx)+kron(diag(ones(1,NH-1),-1),Tx');
hx1=kron(eye(NH),Ty);
hy1=kron(eye(NH),Tz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=length(h01); c=h*N;  L=c;  c1=h*NV; L1=c1;
u1=[1:c1];u2=[c1+1:2*c1];u3=[2*c1+1:c+2*c1];u4=[c+2*c1+1:c+3*c1];u5=[c+3*c1+1:c+4*c1];u6=[c+4*c1+1:2*c+4*c1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%
%左右端口自能；
sigma1=-j*eye(L);  sigma2=-j*eye(L1);
L1L=sigma1;  L1R=sigma1; L1U=sigma2;L1D=sigma2;
parfor i1=1:q
    tic
    for i2=1:time,
        i1,
        E0=Energy(i1);   E=E0+dt*j;
        HE=kron(eye(L),E0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %L1L一层左端口;L1R一层右端口;L1U一层上端口;L1D一层下端口;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %第一三层计算表面格林函数
        HL101=kron(eye(N),hy1);
        HL100=kron(diag(ones(1,N-1),1),hx1)+kron(diag(ones(1,N-1),-1),hx1')+kron(eye(N),h01);
        H01=HL101;
        H00=HL100-HE;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        grii=inv(-H00-L1L);gr1i=grii;gri1=grii;gr11=grii;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k=1:K,
            grii=inv(W*(diag(rand(1,L)-0.5))-H00-H01*grii*H01'); gr1i=gr1i*H01'*grii;
            gr11=gr11+gr1i*H01*gri1; gri1=grii*H01*gri1;
        end
        for k=1:NV,
            grii=inv(W*(diag(rand(1,L)-0.5))-H00-H01*grii*H01'); gr1i=gr1i*H01'*grii;
            gr11=gr11+gr1i*H01*gri1; gri1=grii*H01*gri1;
            gr11=[grii(1:h,1:h),gri1(1:h,:),grii(1:h,h*N-h+1:h*N);gr1i(:,1:h),gr11,gr1i(:,h*N-h+1:h*N);grii(h*N-h+1:h*N,1:h),gri1(h*N-h+1:h*N,:),grii(h*N-h+1:h*N,h*N-h+1:h*N)];
            gr1i=[grii(1:h,:);gr1i;grii(h*N-h+1:h*N,:)];gri1=[grii(:,1:h),gri1,grii(:,h*N-h+1:h*N)];
        end
        for k=1:K,
            grii=inv(W*(diag(rand(1,L)-0.5))-H00-H01*grii*H01'); gr1i=gr1i*H01'*grii;
            gr11=gr11+gr1i*H01*gri1; gri1=grii*H01*gri1;
        end
        for k=1:NV,
            grii=inv(W*(diag(rand(1,L)-0.5))-H00-H01*grii*H01'); gr1i=gr1i*H01'*grii;
            gr11=gr11+gr1i*H01*gri1; gri1=grii*H01*gri1;
            gr11=[grii(1:h,1:h),gri1(1:h,:),grii(1:h,h*N-h+1:h*N);gr1i(:,1:h),gr11,gr1i(:,h*N-h+1:h*N);grii(h*N-h+1:h*N,1:h),gri1(h*N-h+1:h*N,:),grii(h*N-h+1:h*N,h*N-h+1:h*N)];
            gr1i=[grii(1:h,:);gr1i;grii(h*N-h+1:h*N,:)];gri1=[grii(:,1:h),gri1,grii(:,h*N-h+1:h*N)];
        end
        for k=1:K,
            grii=inv(W*(diag(rand(1,L)-0.5))-H00-H01*grii*H01'); gr1i=gr1i*H01'*grii;
            gr11=gr11+gr1i*H01*gri1; gri1=grii*H01*gri1;
        end
        grii=inv(-H00-H01*grii*H01'); gr1i=gr1i*H01'*grii;
        gr11=gr11+gr1i*H01*gri1; gri1=grii*H01*gri1;
        Lgr=[gr11,gr1i;gri1,grii];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        Sigma=zeros(2*c+4*c1);
        Sigma(u1,u1)=L1U;
        Sigma(u2,u2)=L1U;
        Sigma(u4,u4)=L1D;
        Sigma(u5,u5)=L1D;
        Sigma(u6,u6)=L1R;
        LGr1=inv(inv(Lgr)-Sigma);
        
        %格林函数
        GA12=LGr1(u3,u2);GA13=LGr1(u3,u1);GA14=LGr1(u3,u6);GA15=LGr1(u3,u5);GA16=LGr1(u3,u4);
        GA21=LGr1(u2,u3);GA23=LGr1(u2,u1);GA24=LGr1(u2,u6);GA25=LGr1(u2,u5);GA26=LGr1(u2,u4);
        GA31=LGr1(u1,u3);GA32=LGr1(u1,u2);GA34=LGr1(u1,u6);GA35=LGr1(u1,u5);GA36=LGr1(u1,u4);
        GA41=LGr1(u6,u3);GA42=LGr1(u6,u2);GA43=LGr1(u6,u1);GA45=LGr1(u6,u5);GA46=LGr1(u6,u4);
        GA51=LGr1(u5,u3);GA52=LGr1(u5,u2);GA53=LGr1(u5,u1);GA54=LGr1(u5,u6);GA56=LGr1(u5,u4);
        GA61=LGr1(u4,u3);GA62=LGr1(u4,u2);GA63=LGr1(u4,u1);GA64=LGr1(u4,u6);GA65=LGr1(u4,u5);
        
        %GAMMA矩阵
        G1L=j*(L1L-L1L');G1R=j*(L1R-L1R');G1D=j*(L1D-L1D');G1U=j*(L1U-L1U');
        
        
        %无序下的隧穿系数
        T12(i2)=real(trace(G1L*GA12*G1U*GA12'));T13(i2)=real(trace(G1L*GA13*G1U*GA13'));T14(i2)=real(trace(G1L*GA14*G1R*GA14'));T15(i2)=real(trace(G1L*GA15*G1D*GA15'));T16(i2)=real(trace(G1L*GA16*G1D*GA16'));
        T21(i2)=real(trace(G1U*GA21*G1L*GA21'));T23(i2)=real(trace(G1U*GA23*G1U*GA23'));T24(i2)=real(trace(G1U*GA24*G1R*GA24'));T25(i2)=real(trace(G1U*GA25*G1D*GA25'));T26(i2)=real(trace(G1U*GA26*G1D*GA26'));
        T31(i2)=real(trace(G1U*GA31*G1L*GA31'));T32(i2)=real(trace(G1U*GA32*G1U*GA32'));T34(i2)=real(trace(G1U*GA34*G1R*GA34'));T35(i2)=real(trace(G1U*GA35*G1D*GA35'));T36(i2)=real(trace(G1U*GA36*G1D*GA36'));
        T41(i2)=real(trace(G1R*GA41*G1L*GA41'));T42(i2)=real(trace(G1R*GA42*G1U*GA42'));T43(i2)=real(trace(G1R*GA43*G1U*GA43'));T45(i2)=real(trace(G1R*GA45*G1D*GA45'));T46(i2)=real(trace(G1R*GA46*G1D*GA46'));
        T51(i2)=real(trace(G1D*GA51*G1L*GA51'));T52(i2)=real(trace(G1D*GA52*G1U*GA52'));T53(i2)=real(trace(G1D*GA53*G1U*GA53'));T54(i2)=real(trace(G1D*GA54*G1R*GA54'));T56(i2)=real(trace(G1D*GA56*G1D*GA56'));
        T61(i2)=real(trace(G1D*GA61*G1L*GA61'));T62(i2)=real(trace(G1D*GA62*G1U*GA62'));T63(i2)=real(trace(G1D*GA63*G1U*GA63'));T64(i2)=real(trace(G1D*GA64*G1R*GA64'));T65(i2)=real(trace(G1D*GA65*G1D*GA65'));
    end
    toc
    %无序平均后的隧穿系数
    T12=sum(T12)/time;T13=sum(T13)/time;T14=sum(T14)/time;T15=sum(T15)/time;T16=sum(T16)/time;
    T21=sum(T21)/time;T23=sum(T23)/time;T24=sum(T24)/time;T25=sum(T25)/time;T26=sum(T26)/time;
    T31=sum(T31)/time;T32=sum(T32)/time;T34=sum(T34)/time;T35=sum(T35)/time;T36=sum(T36)/time;
    T41=sum(T41)/time;T42=sum(T42)/time;T43=sum(T43)/time;T45=sum(T45)/time;T46=sum(T46)/time;
    T51=sum(T51)/time;T52=sum(T52)/time;T53=sum(T53)/time;T54=sum(T54)/time;T56=sum(T56)/time;
    T61=sum(T61)/time;T62=sum(T62)/time;T63=sum(T63)/time;T64=sum(T64)/time;T65=sum(T65)/time;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T1=T12+T13+T14+T15+T16;
    T2=T21+T23+T24+T25+T26;
    T3=T31+T32+T34+T35+T36;
    T4=T41+T43+T42+T45+T46;
    T5=T51+T53+T54+T52+T56;
    T6=T61+T63+T64+T65+T62;
    V=1;
    D1=(T32*T25+T2*T35)*(T52*T23+T2*T53)+(T2*T3-T32*T23)*(T52*T25-T2*T5);
    D2=(T32*T26+T2*T36)*(T52*T23+T2*T53)+(T2*T3-T32*T23)*(T52*T26+T2*T56);
    D3=(T32*(T21-T24)-T2*(T34-T31))*(T52*T23+T2*T53)+(T2*T3-T32*T23)*(T52*(T21-T24)+T2*(T51-T54));
    D4=((T62*T23+T2*T63)*(T32*T25+T2*T35)+(T2*T3-T32*T23)*(T62*T25+T2*T65));
    D5=((T62*T23+T2*T63)*(T32*T26+T2*T36)+(T2*T3-T32*T23)*(T62*T26-T2*T6));
    D6=((T62*T23+T2*T63)*(T32*(T21-T24)-T2*(T34-T31))+(T2*T3-T32*T23)*(T62*(T21-T24)+T2*(T61-T64)));
    V6=(D4*D3-D1*D6)/(D1*D5-D4*D2);
    V5=(D1*D2*D6-D1*D3*D5)/(D1*D1*D5-D1*D2*D4);
    V3=((T32*T25+T2*T35)*(D1*D2*D6-D1*D3*D5)+(T32*T26+T2*T36)*(D4*D3-D1*D6)*D1+(T32*T21-T32*T24-T2*T34+T2*T31)*D1*(D1*D5-D2*D4))/((T2*T3-T32*T23)*D1*(D1*D5-D2*D4));
    V2=(T23*V3+T25*V5+T26*V6+(T21-T24)*V)/T2;
    I=(T1*V-T12*V2-T13*V3+T14*V-T15*V5-T16*V6);
    
    T=[-T1,T12,T13,T14,T15,T16;
        T21,-T2,T23,T24,T25,T26;
        T31,T32,-T3,T34,T35,T36;
        T41,T42,T43,-T4,T45,T46;
        T51,T52,T53,T54,-T5,T56;
        T61,T62,T63,T64,T65,-T6];
    TT(i1,:,:)=T;
    
    
    
    RH(i1)=(V2-V6)/I;
    RH2(i1)=(V3-V5)/I;
    RL(i1)=(V2-V3)/I;
end

hold on
plot(Energy,TT(:,1,6)+TT(:,2,6)+TT(:,3,6)+TT(:,4,6)+TT(:,5,6),'r');
plot(Energy,1./RH,'b.-');
axis tight






