function one_reactor_opt
%decision variable T
%26 variables
%FCj-x[1-16];Dv-x[17];ri-x[18-19];ki-x[20-21];Ci-x[21-25];T-x[26]
%FCj-  1   2   3   4 
%    A
%    B
%    P    F(i,j)=x(4*(j-1)+i)
%    R             
clear all;clc;
%9 parameters
V=30;
F0=[52.5,50,0,0];
beta=[1,1,0,0];
%Objective funciton
obj=@(x) x(9)+x(10)+x(11)+x(12);
%linear constraints
Aeq=zeros(17,26);
beq=zeros(17,1);
vij=[-1,-1;-1,0;1,-1;0,1];
%MATERIAL BALANCE
%Reactors
for i=1:4
    Aeq(i,(4*(1-1)+i))=1;
    Aeq(i,(4*(2-1)+i))=-1;
    Aeq(i,18:19)=vij(i,:)*V;
end
%Separation
for i=1:4
    Aeq(4+i,(4*(2-1)+i))=1;%Fi,2-Fi,3-Fi,4=0;Fi,3-betai*Fi,2=0
    Aeq(4+i,(4*(3-1)+i))=-1;
    Aeq(4+i,(4*(4-1)+i))=-1;
    Aeq(8+i,(4*(3-1)+i))=1;%Fi,3-betai*Fi,2=0
    Aeq(8+i,(4*(2-1)+i))=-beta(i);
end

%recycle
%for i=1:4
%    y(12+i)=F0(i)+F(i,3)-F(i,1);%x(1)-x(9)=52.5;x(2)-x(10)=50;x(3)-x(11)=0;x(4)-x(12)=0;
%end
for i=1:4
    beq(12+i)=F0(i);%Fi,1-Fi,3=Fi,0
    Aeq(12+i,(4*(1-1)+i))=1;
    Aeq(12+i,(4*(3-1)+i))=-1;
end
%Define Dv
Aeq(17,17)=1;
vm(1)=56.11/620;
vm(2)=58.12/593.4;
vm(3)=114.23/690;
vm(4)=170.33/752;
for i=1:4
    Aeq(17,4*(2-1)+i)=-vm(i);
end
%fmincon input
A=[];b=[];%no inequalities
for i=1:25
    lb(i)=-1e-7;
end
lb(26)=233;
for i=1:25
    ub(i)=1e5;
end
ub(26)=273;
%obtain x0
x0_fsolve=ones(1,25);
[x_fsolve,fval]=fsolve(@CSTR_Structure1,x0_fsolve,[],V);
for i=1:25
    x0(i)=x_fsolve(i);
end
x0(26)=263;
[x,fval]=fmincon(obj,x0,A,b,Aeq,beq,lb,ub,@nonlcon)
F_final=ones(4,5);
for i=1:4
    F_final(i)=F0(i);
end
for i=5:20
    F_final(i)=x(i-4);
end
fprintf('\t0\t\t1\t\t2\t\t3\t\t4\n')
disp(F_final)
fprintf('V=%3.3f\n',V)
fprintf('decision variable T=%3.3f\n',x(26))
fprintf('optimal recycle is %6.5f\n',fval)
function [c,ceq]=nonlcon(x)
%reaction nonlinear constraints
k10=1.66e9*3600;
k20=4.16e12*3600;
R=8.314;
E1=6.5e4;
E2=8.1e4;
ceq(1)=x(22)-x(5)/x(17);
ceq(2)=x(23)-x(6)/x(17);
ceq(3)=x(24)-x(7)/x(17);
ceq(4)=x(25)-x(8)/x(17);
ceq(5)=x(20)-k10*exp(-E1/R/x(26));
ceq(6)=x(21)-k20*exp(-E2/R/x(26));
ceq(7)=x(18)-x(20)*x(22)*x(23);
ceq(8)=x(19)-x(21)*x(22)*x(24);
c=[];

%reactions


%y(18)=k1-k10*exp(-E1/R/T);%x(20)-k10*exp(-E1/R/x(26))
%y(20)=k2-k20*exp(-E2/R/T);%x(21)-k20*exp(-E2/R/x(26))
%y(17)=r1-k1*C(1)*C(2);%x(18)-x(20)*x(22)*x(23)
%y(19)=r2-k2*C(1)*C(3);%x(19)-x(21)*x(22)*x(24)
%all the components are in liquid phase not gaseous
%for i=1:4
%    y(20+i)=C(i)-F(i,2)/Dv;%x(22)-x(5)/x(17);x(23)-x(6)/x(17);x(24)-x(7)/x(17);x(25)-x(8)/x(17);
%end
%y(25)=Dv-VA*F(1,2)-VB*F(2,2)-VP*F(3,2)-VR*F(4,2);%x(17)-VA*x(5)-VB*x(6)-VP*x(7)-VR*x(8)










