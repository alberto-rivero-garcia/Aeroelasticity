%This script will extract the points for a spline curve that will detail
%the separation between the safe and unsafe zones. 

%The spline is based on Ktheta due to the shape of the curve. The logic is
%simple: look for the last Nktheta points in Refined_Grid: these must be
%the most refined points. If they are unsafe, add them to the spline.
%However, if they are safe, add the point to their right to the spline. 


%Load refined grid and compute NKtheta: 
load('Refined_Grid.mat')

NKtheta=length(unique(Refined_Grid(:,1)));
Npoints=length(Refined_Grid(:,1));

%Introduce physical values: 
c=0.537;            %[m]
Ktheta_nom=33032;   %[Nm/rad]

%Extract values of ycg and Ktheta for the spline. 
Spline_Ktheta_v=nan(1,NKtheta);
Spline_Ycg_v=nan(1,NKtheta);

for i=1:NKtheta
    point=Npoints-NKtheta+i;
    if Refined_Grid(point,5) %point is safe
        right_point=Refined_Grid(point,4);
        Spline_Ktheta_v(i)=Refined_Grid(right_point,1);
        Spline_Ycg_v(i)=Refined_Grid(right_point,2);
    else
        Spline_Ktheta_v(i)=Refined_Grid(point,1);
        Spline_Ycg_v(i)=Refined_Grid(point,2);
    end
end

%Prepare vector for spline: 
nplot=1000;
vplot=linspace(0,100,nplot);

%Plot resulting grid and spline: 
figure
clf
axis([24,40,0,100])
hold on
for k=1:Npoints
    if Refined_Grid(k,5)
        scatter(100*Refined_Grid(k,2)/c,100*Refined_Grid(k,1)/Ktheta_nom,'g*')
    else
        scatter(100*Refined_Grid(k,2)/c,100*Refined_Grid(k,1)/Ktheta_nom,'r*')
    end
end
plot(interp1(100*Spline_Ktheta_v/Ktheta_nom,100*Spline_Ycg_v/c,vplot,"cubic"),vplot,'k')
xlabel('Blade CoG [\% Chord]','Interpreter','latex')
ylabel('Control Chain Stiffness [\% Nominal]','Interpreter','latex')

%% Define quantities to be saved for computation of comparison plots: 
Ex3_Spline_Ktheta_v=Spline_Ktheta_v;
Ex3_Spline_Ycg_v=Spline_Ycg_v;

save('Ex3_Spline','Ex3_Spline_Ycg_v',"Ex3_Spline_Ktheta_v")