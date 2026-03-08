%This script computes flutter velocities for a sweep of xmhat values. (0 to
%1). It then plots them and compares them to the baseline case (xmhat=0.33)
%and the 1.15*baseline condition given by exercise 5.d)
addpath("Auxiliary Functions\")

Ne=15;
Nd=6; 
Display=0;
Umax=300;
Nu=2000;
Nx=20;

xmhatv=linspace(0,1,Nx);
UFv=nan(size(xmhatv));

for i=1:Nx
    UFv(i)=PK_TipMassGoland(Ne,Nd,Umax,Nu,xmhatv(i),Display);
end

UFBaseline=PK_TipMassGoland(Ne,Nd,Umax,Nu,0.33,Display);

figure
plot(xmhatv,UFBaseline*ones(size(xmhatv)),'k--','DisplayName','Elastic axis')
hold on 
plot(xmhatv,1.15*UFBaseline*ones(size(xmhatv)),'r--','DisplayName','1.15x Elastic axis')
plot(xmhatv,UFv)
grid on
grid minor
xlabel('$\frac{x_m}{c} \; [-]$','Interpreter','latex')
ylabel('$U_F \; [m/s]$','Interpreter','latex')

legend
