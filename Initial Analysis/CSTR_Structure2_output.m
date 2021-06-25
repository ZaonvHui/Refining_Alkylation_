clear all;clc;
x0=ones(1,26);
%x0=[1,1,0,0,1,1,1,1,1,1,0,0,0,0,1,1,0.1,1,1,1,1,1,1,1,1,1];
[x,fval]=fsolve(@CSTR_Structure2,x0)
F=ones(4,4);
for i=1:16
    F(i)=x(i);
end
F0=[52.5,x(26),0,0];
F_final=ones(4,5);
for i=1:4
    F_final(i)=F0(i);
end
for i=5:20
    F_final(i)=F(i-4);
end

fprintf('\t0\t1\t2\t3\t4\n')
for i=1:4
    for j=1:5
    fprintf('\t%5.3f',F_final(i,j))
    end
        fprintf('\n')
end
if (F_final(1,5)>1e-9) | (F_final(1,5)<-1e-9)
    fprintf('no solution achieved\n')
    conversion=0
else
conversion=(F_final(3,3)+F_final(4,3))/F_final(1,2)
end
