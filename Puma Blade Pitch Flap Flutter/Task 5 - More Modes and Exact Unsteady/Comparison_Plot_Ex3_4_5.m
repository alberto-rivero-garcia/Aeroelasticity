%This script generates a plot which compares the stability maps obtained in
%exercises 3, 4 and 5

%Load spline vectors for exercises 3, 4 and 5
load("Ex3_Spline.mat")
load("Ex4_Spline.mat")
load("Ex5_Spline.mat")

%Introduce physical values: 
c=0.537;            %[m]
Ktheta_nom=33032;   %[Nm/rad]

%Prepare vector for spline plots: 
nplot=1000;
vplot=linspace(0,100,nplot);

%Plot resulting grid and spline: 
figure
clf
axis([24,40,0,100])
hold on
plot(interp1(100*Ex3_Spline_Ktheta_v/Ktheta_nom,100*Ex3_Spline_Ycg_v/c,vplot,"cubic"),vplot,'DisplayName','Two Modes, Reference Section')
plot(interp1(100*Ex4_Spline_Ktheta_v/Ktheta_nom,100*Ex4_Spline_Ycg_v/c,vplot,"cubic"),vplot,'DisplayName','Six Modes, Reference Section')
plot(interp1(100*Ex5_Spline_Ktheta_v/Ktheta_nom,100*Ex5_Spline_Ycg_v/c,vplot,"cubic"),vplot,'DisplayName','Six Modes, Exact Unsteady')
xlabel('Blade CoG [\% Chord]','Interpreter','latex')
ylabel('Control Chain Stiffness [\% Nominal]','Interpreter','latex')
legend()