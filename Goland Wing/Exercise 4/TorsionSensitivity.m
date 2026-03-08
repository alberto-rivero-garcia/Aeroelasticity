%This script performs a sensitivity analysis of the flutter velocity wrt
%the structural torsion stiffness: 
addpath("Auxiliary Functions\")
%Setup basic values: (I'll fine tune them later on)
Ne=15;
Nd=6; 
Ns=5; %number of points for sensitivity analysis
Umax=200;
Nu=1000;

%Specify parameters: 
params.EI=9.77*1e6;    %[N*m2]
params.ghat=0.43;      %[-]
params.Clalpha=2*pi;   %[-]
%Specify nominal GJ
GJnom=987600;      %[N m^2]
GJv=linspace(0.98*GJnom,1.02*GJnom,Ns);

%Prepare iteration: 
UFv=nan(size(GJv));

for i=1:Ns
    params.GJ=GJv(i);
    UFv(i)=SensitivityPKGoland(Ne,Nd,Umax,Nu,params);
end

params.GJ=GJnom;
UFnom=SensitivityPKGoland(Ne,Nd,Umax,Nu,params);

%Plot results: 
GJpercent=linspace(-2,2,Ns);
UFpercent=(UFv-UFnom)/UFnom*100;

figure(1)
plot(GJpercent,UFpercent)
hold on
grid minor
xlabel("$\Delta GJ \;[\%]$",'Interpreter','latex')
ylabel("$\Delta U_F \;[\%]$",'Interpreter','latex')

%REQUIRES IMPROVEMENT OF UF INTERPOLATION TECHNIQUE IN PK METHOD -> APPLY
%IT EVERYWHERE. 
