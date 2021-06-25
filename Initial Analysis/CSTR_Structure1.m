function y=CSTR_Structure1(x,V)
%25 equations & 25 variables
T=263;
F0=[52.5,50,0,0];
beta=[1,1,0,0];
F=ones(4,4);
for i=1:16
    F(i)=x(i);
end
Dv=x(17);
r1=x(18);r2=x(19);
k1=x(20);k2=x(21);
for i=1:4
    C(i)=x(21+i);
end
%CA=x(22);CB=x(23);CP=x(24);CR=x(25)
%mass balance
%A
y(1)=F(1,1)-F(1,2)+V*(-r1-r2);
%B
y(2)=F(2,1)-F(2,2)+V*(-r1);
%P
y(3)=F(3,1)-F(3,2)+V*(r1-r2);
%R
y(4)=F(4,1)-F(4,2)+V*r2;
%separation point
for i=1:4
    y(4+i)=F(i,2)-F(i,3)-F(i,4);
    y(8+i)=F(i,3)-beta(i)*F(i,2);
end
%recycle
for i=1:4
    y(12+i)=F0(i)+F(i,3)-F(i,1);
end
%reactions
k10=1.66e9*3600;
k20=4.16e12*3600;
R=8.314;
E1=6.5e4;
E2=8.1e4;
y(18)=k1-k10*exp(-E1/R/T);
y(20)=k2-k20*exp(-E2/R/T);
y(17)=r1-k1*C(1)*C(2);
y(19)=r2-k2*C(1)*C(3);
%all the components are in liquid phase not gaseous
VA=56.11/620;
VB=58.12/593.4;
VP=114.23/690;
VR=170.33/752;
for i=1:4
    y(20+i)=C(i)-F(i,2)/Dv;
end
y(25)=Dv-VA*F(1,2)-VB*F(2,2)-VP*F(3,2)-VR*F(4,2);

