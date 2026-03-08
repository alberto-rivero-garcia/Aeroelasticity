%This script performs a sensitivity analysis of the flutter velocity wrt
%the lift curve slope:  
addpath("Auxiliary Functions\")
%Setup basic values: (I'll fine tune them later on)
Ne=15;
Nd=6; 
Ns=5; %number of points for sensitivity analysis
Umax=200;
Nu=1000;

%Specify parameters: 
params.EI=9.77*1e6;    %[N*m2]
params.GJ=987600;      %[N m^2]
params.ghat=0.43;      %[-]

%Specify nominal Clalpha:
Clalphanom=2*pi;
Clalphav=linspace(0.98*Clalphanom,1.02*Clalphanom,Ns);

%Prepare iteration: 
UFv=nan(size(Clalphav));

for i=1:Ns
    params.Clalpha=Clalphav(i);
    UFv(i)=SensitivityPKGoland(Ne,Nd,Umax,Nu,params);
end

params.Clalpha=Clalphanom;
UFnom=SensitivityPKGoland(Ne,Nd,Umax,Nu,params);

%Plot results: 
Clalphapercent=linspace(-2,2,Ns);
UFpercent=(UFv-UFnom)/UFnom*100;

figure(1)
plot(Clalphapercent,UFpercent)
hold on
grid minor
xlabel("$\Delta C_{l \alpha} \;[\%]$",'Interpreter','latex')
ylabel("$\Delta U_F \;[\%]$",'Interpreter','latex')

%REQUIRES IMPROVEMENT OF UF INTERPOLATION TECHNIQUE IN PK METHOD -> APPLY
%IT EVERYWHERE. 
