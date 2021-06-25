function three_reactor_C4_opt
%decision variable s1 s2 s3 V1 V2 V3 T
%90 variables
%FCj-x[1-54];Dv-x[55-57];ri-x[58-63];ki-x[64-65];Ci-x[66-83];T-x[84];Vi-x[85-87];si-x[88-90];
%FCj-  1 2 3 4 0a1(5) 0a2(6) 0a3(7) 20(8) 21(9)
%    A
%    B
%    P
%    R              F(i,j)=x(6*(j-1)+i)
%    I
%    N
%ri-x[58-63]r11 r12 r21 r22 r31 r32;
%Ci-x[66-83]CA1 CB1 CP1 CR1 CI1 CN1 CA2 CB2 CP2 CR2 CI2 CN2 CA3 CB3 CP3 CR3 CI3 CN3;clear all;clc;
clear all;clc;
%18 parameters
beta=[1,1,0,0,0,0];
F0a=[52.5,20,0,0,5.5,27];
F0b=[0,30,0,0,3,27];
%Objective funciton
obj=@(x) x(13)+x(14)+x(15)+x(16)+x(17)+x(18);
%linear constraints
Aeq=zeros(40,90);
beq=zeros(40,1);
%MATERIAL BALANCE
%Separation
beta=[1,1,0,0,0,0];
F0a=[52.5,20,0,0,5.5,27];
F0b=[0,30,0,0,3,27];
for i=1:6
    Aeq(i,(6*(2-1)+i))=1;%Fi,2-Fi,3-Fi,4=0;Fi,3-betai*Fi,2=0
    Aeq(i,(6*(3-1)+i))=-1;
    Aeq(i,(6*(4-1)+i))=-1;
    Aeq(6+i,(6*(3-1)+i))=1;
    Aeq(6+i,(6*(2-1)+i))=-beta(i);
    
    
    Aeq(12+i,(6*(5-1)+i))=1;%Fi,0a1-Fi,0a*s1=0;Fi,0a2-Fi,0a*s2=0;Fi,0a3-Fi,0a*s3=0
    Aeq(12+i,88)=-F0a(i);
    Aeq(18+i,(6*(6-1)+i))=1;
    Aeq(18+i,89)=-F0a(i);
    Aeq(24+i,(6*(7-1)+i))=1;
    Aeq(24+i,90)=-F0a(i);
end

%Recycle
for i=1:6
    beq(30+i)=F0b(i);%Fi,1-Fi,3=Fi,0b
    Aeq(30+i,(6*(1-1)+i))=1;
    Aeq(30+i,(6*(3-1)+i))=-1;
end
%Define Dv
Aeq(37,55)=1;
Aeq(38,56)=1;
Aeq(39,57)=1;
vm(1)=56.11/620;
vm(2)=58.12/593.4;
vm(3)=114.23/690;
vm(4)=170.33/752;
vm(5)=44.1/493;
vm(6)=58.12/573;
for i=1:6
    Aeq(37,6*(8-1)+i)=-vm(i);%Dv(1)-sum(Fi,20*vm(i))=0
    Aeq(38,6*(9-1)+i)=-vm(i);%Dv(2)-sum(Fi,21*vm(i))=0
    Aeq(39,6*(2-1)+i)=-vm(i);%Dv(3)-sum(Fi,2*vm(i))=0
end
Aeq(40,88)=1;Aeq(40,89)=1;Aeq(40,90)=1;beq(40)=1;

%fmincon input
A=[];b=[];
for i=1:90
    lb(i)=-1e-7;
end
lb(84)=233;lb(85:90)=0;
for i=1:87
    ub(i)=1e5;
end
ub(84)=273;ub(85:87)=100;ub(88:90)=1;
%obtain x0
x0_fsolve=ones(1,84);%?
[x_fsolve,fval]=fsolve(@three_reactor_output,x0_fsolve);
for i=1:83
    x0(i)=x_fsolve(i);
end
x0(84)=263;
x0(85:87)=20;
x0(88:90)=1/3;
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
fprintf('decision variable s1=%3.3f\n',x(88))
fprintf('decision variables2=%3.3f\n',x(89))
fprintf('decision variable s3=%3.3f\n',x(90))
fprintf('decision variable V1=%3.3f\n',x(85))
fprintf('decision variable V2=%3.3f\n',x(86))
fprintf('decision variable V3=%3.3f\n',x(87))
fprintf('decision variable T=%3.3f\n',x(84))
fprintf('optimal recycle is %6.5f\n',fval)
function [c,ceq]=nonlcon(x)
c=[];
k10=1.66e9*3600;
k20=4.16e12*3600;
R=8.314;
E1=6.5e4;
E2=8.1e4;
ceq(1)=x(64)-k10*exp(-E1/R/x(84));
ceq(2)=x(65)-k20*exp(-E2/R/x(84));
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

ceq(27)=x(6*(1-1)+1)+x(6*(5-1)+1)-x(6*(8-1)+1)+x(85)*(-x(58)-x(59));
ceq(28)=x(6*(1-1)+2)+x(6*(5-1)+2)-x(6*(8-1)+2)+x(85)*(-x(58));
ceq(29)=x(6*(1-1)+3)+x(6*(5-1)+3)-x(6*(8-1)+3)+x(85)*(x(58)-x(59));
ceq(30)=x(6*(1-1)+4)+x(6*(5-1)+4)-x(6*(8-1)+4)+x(85)*(x(59));
ceq(31)=x(6*(1-1)+5)+x(6*(5-1)+5)-x(6*(8-1)+5);
ceq(32)=x(6*(1-1)+6)+x(6*(5-1)+6)-x(6*(8-1)+6);
ceq(33)=x(6*(8-1)+1)+x(6*(6-1)+1)-x(6*(9-1)+1)+x(86)*(-x(60)-x(61));
ceq(34)=x(6*(8-1)+2)+x(6*(6-1)+2)-x(6*(9-1)+2)+x(86)*(-x(60));
ceq(35)=x(6*(8-1)+3)+x(6*(6-1)+3)-x(6*(9-1)+3)+x(86)*(x(60)-x(61));
ceq(36)=x(6*(8-1)+4)+x(6*(6-1)+4)-x(6*(9-1)+4)+x(86)*(x(61));
ceq(37)=x(6*(8-1)+5)+x(6*(6-1)+5)-x(6*(9-1)+5);
ceq(38)=x(6*(8-1)+6)+x(6*(6-1)+6)-x(6*(9-1)+6);
ceq(39)=x(6*(9-1)+1)+x(6*(7-1)+1)-x(6*(2-1)+1)+x(87)*(-x(62)-x(63));
ceq(40)=x(6*(9-1)+2)+x(6*(7-1)+2)-x(6*(2-1)+2)+x(87)*(-x(62));
ceq(41)=x(6*(9-1)+3)+x(6*(7-1)+3)-x(6*(2-1)+3)+x(87)*(x(62)-x(63));
ceq(42)=x(6*(9-1)+4)+x(6*(7-1)+4)-x(6*(2-1)+4)+x(87)*(x(63));
ceq(43)=x(6*(9-1)+5)+x(6*(7-1)+5)-x(6*(2-1)+5);
ceq(44)=x(6*(9-1)+6)+x(6*(7-1)+6)-x(6*(2-1)+6);


%FCj-x[1-54];Dv-x[55-57];ri-x[58-63]r11 r12 r21 r22 r31 r32;ki-x[64-65];Ci-x[66-83]CA1 CB1 CP1;si-x[84];T-x[85];s1 s2-x[86 87];Vi-x[88-90];

