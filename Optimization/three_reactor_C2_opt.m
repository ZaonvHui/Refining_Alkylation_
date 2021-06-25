function three_reactor_C2_opt
%Equal volume-decision variable s1 s2 T
%87 variables
%FCj-x[1-54];Dv-x[55-57];ri-x[58-63];ki-x[64-65];Ci-x[66-83];s3-x[84];T-x[85];s1s2-x[86 87]
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
%21 parameters
V1=20;V2=20;V3=20;
beta=[1,1,0,0,0,0];
F0a=[52.5,20,0,0,5.5,27];
F0b=[0,30,0,0,3,27];
%Objective funciton
obj=@(x) x(13)+x(14)+x(15)+x(16)+x(17)+x(18);
%linear constraints
Aeq=zeros(58,87);
beq=zeros(58,1);
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
    Aeq(30+i,86)=-F0a(i);
    Aeq(36+i,(6*(6-1)+i))=1;
    Aeq(36+i,87)=-F0a(i);
    Aeq(42+i,(6*(7-1)+i))=1;
    Aeq(42+i,84)=-F0a(i);
end
Aeq(49,84)=1;Aeq(49,86)=1;Aeq(49,87)=1;beq(49)=1;%s1+s2+s3=1

%Recycle
for i=1:6
    beq(49+i)=F0b(i);%Fi,1-Fi,3=Fi,0b
    Aeq(49+i,(6*(1-1)+i))=1;
    Aeq(49+i,(6*(3-1)+i))=-1;
end
%Define Dv
Aeq(56,55)=1;
Aeq(57,56)=1;
Aeq(58,57)=1;
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

%fmincon input
A=[];b=[];
for i=1:87
    lb(i)=-1e-7;
end
lb(85)=233;lb(84)=0;lb(86:87)=0;
for i=1:87
    ub(i)=1e5;
end
ub(85)=273;ub(84)=1;ub(86:87)=1;
%obtain x0
x0_fsolve=ones(1,84);%?
[x_fsolve,fval]=fsolve(@three_reactor_output,x0_fsolve);
for i=1:83
    x0(i)=x_fsolve(i);
end
x0(85)=263;
x0(84)=1/3;x0(86:87)=1/3;
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
fprintf('V1=%3.3f\n',V1)
fprintf('V2=%3.3f\n',V2)
fprintf('V3=%3.3f\n',V3)
fprintf('decision variable s1=%3.3f\n',x(86))
fprintf('decision variable s2=%3.3f\n',x(87))
fprintf('s3=%3.3f\n',x(84))
fprintf('decision variable T=%3.3f\n',x(85))
fprintf('optimal recycle is %6.5f\n',fval)
function [c,ceq]=nonlcon(x)
c=[];
k10=1.66e9*3600;
k20=4.16e12*3600;
R=8.314;
E1=6.5e4;
E2=8.1e4;
ceq(1)=x(64)-k10*exp(-E1/R/x(85));
ceq(2)=x(65)-k20*exp(-E2/R/x(85));
for i=1:6
    ceq(2+i)=x(65+i)-x(6*(8-1)+i)/x(55);
    ceq(8+i)=x(71+i)-x(6*(9-1)+i)/x(56);
    ceq(14+i)=x(77+i)-x(6*(2-1)+i)/x(57);
end
ceq(21)=x(58)-x(64)*x(66)*x(67);
ceq(22)=x(59)-x(65)*x(66)*x(68);
ceq(23)=x(60)-x(64)*x(72)*x(73);
ceq(24)=x(61)-x(65)*x(72)*x(74);
ceq(25)=x(62)-x(64)*x(78)*x(79);
ceq(26)=x(63)-x(65)*x(78)*x(80);

