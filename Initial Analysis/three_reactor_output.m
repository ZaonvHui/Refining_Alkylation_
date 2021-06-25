function three_reactor_output
clear all;clc;
x0=ones(1,84);%?
[x,fval]=fsolve(@three_reactor,x0);
disp(fval)
F0a=[52.5,20,0,0,5.5,27];
F0b=[0,30,0,0,3,27];
s1=1/3;s2=1/3;V1=20;V2=20;V3=20;T=263;
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
fprintf('s1 is %3.3f\n',s1)
fprintf('s2 is %3.3f\n',s2)
fprintf('s3 is %3.3f\n',x(84))
fprintf('V1 is %3.3f\n',V1)
fprintf('V2 is %3.3f\n',V2)
fprintf('V3 is %3.3f\n',V3)
fprintf('t is %3.3f\n',T)

recycle=sum(F_final(:,3))
