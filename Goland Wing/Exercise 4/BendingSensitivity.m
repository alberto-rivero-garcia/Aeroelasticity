%This script performs a sensitivity analysis of the flutter velocity wrt
%the structural bending stiffness: 
addpath("Auxiliary Functions\")
%Setup basic values: (I'll fine tune them later on)
Ne=15;
Nd=6; 
Ns=5; %number of points for sensitivity analysis
Umax=200;
Nu=1000;

%Specify parameters: 
params.GJ=987600;      %[N m^2]
params.ghat=0.43;      %[-]
params.Clalpha=2*pi;   %[-]
%Specify nominal EI
EInom=9.77*1e6;    %[N*m2]
EIv=linspace(0.98*EInom,1.02*EInom,Ns);

%Prepare iteration: 
UFv=nan(size(EIv));

for i=1:Ns
    params.EI=EIv(i);
    UFv(i)=SensitivityPKGoland(Ne,Nd,Umax,Nu,params);
end

params.EI=EInom;
UFnom=SensitivityPKGoland(Ne,Nd,Umax,Nu,params);

%Plot results: 
EIpercent=linspace(-2,2,Ns);
UFpercent=(UFv-UFnom)/UFnom*100;

figure(1)
plot(EIpercent,UFpercent)
hold on
grid minor
xlabel("$\Delta EI \;[\%]$",'Interpreter','latex')
ylabel("$\Delta U_F \;[\%]$",'Interpreter','latex')

%REQUIRES IMPROVEMENT OF UF INTERPOLATION TECHNIQUE IN PK METHOD -> APPLY
%IT EVERYWHERE. 
