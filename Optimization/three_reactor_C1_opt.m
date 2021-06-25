function three_reactor_C1_opt
%Equal volume & equal split-decision variable T
%84 variables
%FCj-x[1-54];Dv-x[55-57];ri-x[58-63];ki-x[64-65];Ci-x[66-83];T-x[84].
%FCj-  1 2 3 4 0a1(5) 0a2(6) 0a3(7) 20(8) 21(9)
%    A
%    B
%    P
%    R              F(i,j)=x(6*(j-1)+i)
%    I
%    N
%ri-x[58-63]r11 r12 r21 r22 r31 r32;
%Ci-x[66-83]CA1 CB1 CP1 CR1 CI1 CN1 CA2 CB2 CP2 CR2 CI2 CN2 CA3 CB3 CP3 CR3 CI3 CN3;
clear all;clc;
%24 parameters
s1=1/3;s2=1/3;s3=1/3;V1=20;V2=20;V3=20;
beta=[1,1,0,0,0,0];
F0a=[52.5,20,0,0,5.5,27];
F0b=[0,30,0,0,3,27];
%Objective funciton
obj=@(x) x(13)+x(14)+x(15)+x(16)+x(17)+x(18);
%linear constraints
Aeq=zeros(57,84);
beq=zeros(57,1);
vij=[-1,-1;-1,0;1,-1;0,1;0,0;0,0];
%MATERIAL BALANCE
%Reactors
for i=1:6
    Aeq(i,(6*(1-1)+i))=1;
    Aeq(i,(6*(5-1)+i))=1;
    Aeq(i,(6*(8-1)+i))=-1;
    Aeq(i,58:59)=vij(i,:)*V1;
end
for i=1:6
    Aeq(6+i,(6*(8-1)+i))=1;
    Aeq(6+i,(6*(6-1)+i))=1;
    Aeq(6+i,(6*(9-1)+i))=-1;
    Aeq(6+i,60:61)=vij(i,:)*V2;
end
for i=1:6
    Aeq(12+i,(6*(9-1)+i))=1;
    Aeq(12+i,(6*(7-1)+i))=1;
    Aeq(12+i,(6*(2-1)+i))=-1;
    Aeq(12+i,62:63)=vij(i,:)*V3;
end
%Separation
for i=1:6
    Aeq(18+i,(6*(2-1)+i))=1;%Fi,2-Fi,3-Fi,4=0;Fi,3-betai*Fi,2=0
    Aeq(18+i,(6*(3-1)+i))=-1;
    Aeq(18+i,(6*(4-1)+i))=-1;
    Aeq(24+i,(6*(3-1)+i))=1;
    Aeq(24+i,(6*(2-1)+i))=-beta(i);
    
    Aeq(30+i,(6*(5-1)+i))=1;%Fi,0a1=Fi,0a*s1;Fi,0a2=Fi,0a*s2;Fi,0a3=Fi,0a*s3
    beq(30+i)=F0a(i)*s1;
    Aeq(36+i,(6*(6-1)+i))=1;
    beq(36+i)=F0a(i)*s2;
    Aeq(42+i,(6*(7-1)+i))=1;
    beq(42+i)=F0a(i)*s3;
end
%Recycle
for i=1:6
    beq(48+i)=F0b(i);%Fi,1-Fi,3=Fi,0b
    Aeq(48+i,(6*(1-1)+i))=1;
    Aeq(48+i,(6*(3-1)+i))=-1;
end
%Define Dv
Aeq(55,55)=1;
Aeq(56,56)=1;
Aeq(57,57)=1;
vm(1)=56.11/620;
vm(2)=58.12/593.4;
vm(3)=114.23/690;
vm(4)=170.33/752;
vm(5)=44.1/493;
vm(6)=58.12/573;
for i=1:6
    Aeq(55,6*(8-1)+i)=-vm(i);%Dv(1)-sum(Fi,20*vm(i))=0
    Aeq(56,6*(9-1)+i)=-vm(i);%Dv(2)-sum(Fi,21*vm(i))=0
    Aeq(57,6*(2-1)+i)=-vm(i);%Dv(3)-sum(Fi,2*vm(i))=0
end
%fmincon input
A=[];b=[];%no inequalities
for i=1:84
    lb(i)=0;
end
lb(84)=233;
for i=1:84
    ub(i)=1e5;
end
ub(84)=273;
%obtain x0
x0_fsolve=ones(1,84);%?
[x_fsolve,fval]=fsolve(@three_reactor_output,x0_fsolve);
for i=1:83
    x0(i)=x_fsolve(i);
end
x0(84)=263;
[x,fval]=fmincon(obj,x0,A,b,Aeq,beq,lb,ub,@nonlcon);
F_final=ones(6,11);
for i=1:6
    F_final(i)=F0a(i);
end
for i=7:12
    F_final(i)=F0b(i-6);
end
for i=13:54
    F_final(i)=x(i-12);
end
fprintf('\toa\t\tob\t\t1\t\t2\t\t3\t\t4\t\t0a1\t\t0a2\t\t0a3\t\t20\t\t21\n')
disp(F_final)
fprintf('s1=%3.3f\n',s1)
fprintf('s2=%3.3f\n',s2)
fprintf('s3=%3.3f\n',s2)
fprintf('V1=%3.3f\n',V1)
fprintf('V2=%3.3f\n',V2)
fprintf('V3=%3.3f\n',V3)
fprintf('decision variable T=%3.3f\n',x(84))
fprintf('optimal recycle is %6.5f\n',fval)

%nonlinear constraints
function [c,ceq]=nonlcon(x)
c=[];
k10=1.66e9*3600;
k20=4.16e12*3600;
R=8.314;
E1=6.5e4;
E2=8.1e4;
ceq(1)=x(64)-k10*exp(-E1/R/x(84));%ki-ki0*exp(-Ei/R/T);
ceq(2)=x(65)-k20*exp(-E2/R/x(84));
%Cij
for i=1:6
    ceq(2+i)=x(65+i)-x(6*(8-1)+i)/x(55);
    ceq(8+i)=x(71+i)-x(6*(9-1)+i)/x(56);
    ceq(14+i)=x(77+i)-x(6*(2-1)+i)/x(57);
end
%Rate
ceq(21)=x(58)-x(64)*x(66)*x(67);
ceq(22)=x(59)-x(65)*x(66)*x(68);
ceq(23)=x(60)-x(64)*x(72)*x(73);
ceq(24)=x(61)-x(65)*x(72)*x(74);
ceq(25)=x(62)-x(64)*x(78)*x(79);
ceq(26)=x(63)-x(65)*x(78)*x(80);

