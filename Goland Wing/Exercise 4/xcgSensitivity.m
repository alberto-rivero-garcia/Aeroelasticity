%This script performs a sensitivity analysis of the flutter velocity wrt
%the chordwise position of the center of mass: 
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
params.Clalpha=2*pi;   %[-]

%Specify nominal ghat
c=1.829;        %[m]
xcgnom=0.43*c;  %[m]
ghatnom=0.43;
ghatv=1/c*linspace(0.98*xcgnom,1.02*xcgnom,Ns);

%Prepare iteration: 
UFv=nan(size(ghatv));

for i=1:Ns
    params.ghat=ghatv(i);
    UFv(i)=SensitivityPKGoland(Ne,Nd,Umax,Nu,params);
end

params.ghat=ghatnom;
UFnom=SensitivityPKGoland(Ne,Nd,Umax,Nu,params);

%Plot results: 
ghatpercent=linspace(-2,2,Ns);
UFpercent=(UFv-UFnom)/UFnom*100;

figure(1)
plot(ghatpercent,UFpercent)
hold on
grid minor
xlabel("$\Delta \frac{x_{cg}}{c} \;[\%]$",'Interpreter','latex')
ylabel("$\Delta U_F \;[\%]$",'Interpreter','latex')

%REQUIRES IMPROVEMENT OF UF INTERPOLATION TECHNIQUE IN PK METHOD -> APPLY
%IT EVERYWHERE. 
