function y=three_reactor(x)
%84 equations & 108 variables

%24 parameters
s1=1/3;s2=1/3;V1=20;V2=20;V3=20;T=263;
beta=[1,1,0,0,0,0];
F0a=[52.5,20,0,0,5.5,27];
F0b=[0,30,0,0,3,27];

%84 variables
%FCj-x[1-54];Dv-x[55-57];ri-x[58-63];ki-x[64-65];Ci-x[66-83];s3-x[84].
%FCj-  1 2 3 4 0a1(5) 0a2(6) 0a3(7) 20(8) 21(9)
%    A
%    B
%    P
%    R              F(i,j)=x(6*(j-1)+i)
%    I
%    N
%ri-x[58-63]r11 r12 r21 r22 r31 r32;
%Ci-x[66-83]CA1 CB1 CP1 CR1 CI1 CN1 CA2 CB2 CP2 CR2 CI2 CN2 CA3 CB3 CP3 CR3 CI3 CN3;

%84 equations
Aeq=zeros(58,84);%58 linear equations & 16 nonlinear equations
beq=zeros(58,1);
vij=[-1,-1;-1,0;1,-1;0,1;0,0;0,0];%stoitiometric coefficients
%MATERIAL BALANCE
%Reactor 1 R1:FA,1+FA,5-FA,8+V1*(-r11-r12)=0
for i=1:6
    Aeq(i,(6*(1-1)+i))=1;
    Aeq(i,(6*(5-1)+i))=1;
    Aeq(i,(6*(8-1)+i))=-1;
    Aeq(i,58:59)=vij(i,:)*V1;
end
%Reactor 2 R2:FA,20+FA,0a2-FA,21+V2*(-r21-r22)
for i=1:6
    Aeq(6+i,(6*(8-1)+i))=1;
    Aeq(6+i,(6*(6-1)+i))=1;
    Aeq(6+i,(6*(9-1)+i))=-1;
    Aeq(6+i,60:61)=vij(i,:)*V2;
end
%Reactor 3 R3:FA,21+FA,0a3-FA,2+V3*(-r31-r32)
for i=1:6
    Aeq(12+i,(6*(9-1)+i))=1;
    Aeq(12+i,(6*(7-1)+i))=1;
    Aeq(12+i,(6*(2-1)+i))=-1;
    Aeq(12+i,62:63)=vij(i,:)*V3;
end
%Separation_output&input
for i=1:6
    Aeq(18+i,(6*(2-1)+i))=1;%Fi,2-Fi,3-Fi,4=0;Fi,3-betai*Fi,2=0
    Aeq(18+i,(6*(3-1)+i))=-1;
    Aeq(18+i,(6*(4-1)+i))=-1;
    Aeq(24+i,(6*(3-1)+i))=1;
    Aeq(24+i,(6*(2-1)+i))=-beta(i);
    
    Aeq(30+i,(6*(5-1)+i))=1;%Fi,0a1=Fi,0a*s1;Fi,0a2=Fi,0a*s2;Fi,0a3-Fi,0a*s3=0
    beq(30+i)=F0a(i)*s1;
    Aeq(36+i,(6*(6-1)+i))=1;
    beq(36+i)=F0a(i)*s2;
    Aeq(42+i,(6*(7-1)+i))=1;
    Aeq(42+i,84)=-F0a(i);
end
Aeq(49,84)=1;beq(49)=1-s1-s2;%s3=1-s1-s2

%Recycle
for i=1:6
    beq(49+i)=F0b(i);%Fi,1-Fi,3=Fi,0b
    Aeq(49+i,(6*(1-1)+i))=1;
    Aeq(49+i,(6*(3-1)+i))=-1;
end
%Define Dv
Aeq(56,55)=1;%Dv(1)
Aeq(57,56)=1;%Dv(2)
Aeq(58,57)=1;%Dv(3)
vm(1)=56.11/620;
vm(2)=58.12/593.4;
vm(3)=114.23/690;
vm(4)=170.33/752;
vm(5)=44.1/493;
vm(6)=58.12/573;
for i=1:6
Aeq(56,6*(8-1)+i)=-vm(i);%Dv(1)-sum(Fi,20*vm(i))=0
Aeq(57,6*(9-1)+i)=-vm(i);%Dv(2)-sum(Fi,21*vm(i))=0
Aeq(58,6*(2-1)+i)=-vm(i);%Dv(3)-sum(Fi,2*vm(i))=0
end
%Convert matrix Aeq*x=beq into fsolve equations:y(i)=f(x)=0--Aeq*x-beq=0
for i=1:58
    y(i)=0;
    for j=1:84
        y(i)=y(i)+Aeq(i,j)*x(j);
    end
    y(i)=y(i)-beq(i);
end

%REACTIONS
k10=1.66e9*3600;
k20=4.16e12*3600;
R=8.314;
E1=6.5e4;
E2=8.1e4;
%ki  ki-ki0*exp(-Ei/R/T);
y(59)=x(64)-k10*exp(-E1/R/T);
y(60)=x(65)-k20*exp(-E2/R/T);
%Cij-x[66-83]CA1 CB1 CP1 CR1 CI1 CN1 CA2 CB2 CP2 CR2 CI2 CN2 CA3 CB3 CP3 CR3 CI3 CN3;
for i=1:6
    y(60+i)=x(65+i)-x(6*(8-1)+i)/x(55);%R1
    y(66+i)=x(71+i)-x(6*(9-1)+i)/x(56);%R2
    y(72+i)=x(77+i)-x(6*(2-1)+i)/x(57);%R3
end
%Rate ri-x[58-63]r11 r12 r21 r22 r31 r32;
y(79)=x(58)-x(64)*x(66)*x(67);%r11=k1*CA1*CB1
y(80)=x(59)-x(65)*x(66)*x(68);%r12=k2*CA1*CP1
y(81)=x(60)-x(64)*x(72)*x(73);%r21=k1*CA2*CB2
y(82)=x(61)-x(65)*x(72)*x(74);%r22=k2*CA2*CP2
y(83)=x(62)-x(64)*x(78)*x(79);%r31=k1*CA3*CB3
y(84)=x(63)-x(65)*x(78)*x(80);%r32=k2*CA3*CP3


  




