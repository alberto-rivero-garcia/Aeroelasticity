%This script computes the baseline grid for the stability map. It begins
%with an equispaced 25 element grid and will do the first stability plot of
%task 3. Other routines will then start from this grid and refine it in
%such a way to better detect the boundaries of the safe region. 



%Basic file management
addpath("Data\","FEM\")

tic
%Choose grid points:
Ncg=4; 
NKtheta=10; 

%Define size of Galerkin discretizations and set continuation algorithm
%parameters. 
Ne=20;
Nd=2;
solparams.Nomega=10;
solparams.maxiter=50; %ample to allow for complete convergence everywhere. 
solparams.tol=1e-6;

%Construct discretization: 
c=0.537;            %[m]
Ktheta_nom=33032;   %[Nm/rad]


Ycg_hat_v=linspace(0.24,0.4,Ncg)*c;
Ktheta_v=linspace(0,1,NKtheta)*Ktheta_nom;

%Create point list. Make it so the lowest row is first, then the second and
%so on. Note that Ktheta will be on the vertical and ycg on the horizontal

Baseline_Grid=nan(NKtheta*Ncg,6);

%Grid 1: Ktheta, Grid 2: Ycg_hat, Grid 3: Left Neighbour, Grid 4: Right
%neighbour, Grid 5: safe, Grid 6: convergence

for i=1:NKtheta
    %Points: 
    points=Ncg*(i-1)+(1:Ncg);
    %Update Ktheta
    Baseline_Grid(points,1)=Ktheta_v(i);
    %Update Ycg_hat
    Baseline_Grid(points,2)=Ycg_hat_v';
    %Update left neighbours:
    Baseline_Grid(points(2:end),3)=points(2:end)-1;
    %Update right neighbours: 
    Baseline_Grid(points(1:end-1),4)=points(1:end-1)+1;
end

%Evaluate points on the list and construct a primitive stability map

%Prepare computation of epsilon:  
xe= 0.289;      %[m]
xco=2.080;      %[m]
R=7.490;        %[m]
[~,waypoints_1]=m_Puma([xe,R]);
[~,waypoints_2]=ycg_Puma([xe,R]);
waypoints=sort(unique([waypoints_1;waypoints_2]));
Blade_mass=integral(@(x)m_Puma(x),xe,R,'Waypoints',waypoints_1);
Blade_co_mass=integral(@(x)m_Puma(x),xco,R,'Waypoints',waypoints_1);
Ycg_baseline=integral(@(x)m_Puma(x).*ycg_Puma(x),xe,R,'Waypoints',waypoints)/Blade_mass;

%Compute safety at each point
for i=1:(NKtheta*Ncg)
    %Compute epsilon:
    Ycg_prime=c/4-Baseline_Grid(i,2);
    epsilon=(Ycg_prime-Ycg_baseline)*Blade_mass/Blade_co_mass;
    %Prepare cg input
    f_eps=@(x)epsilon*(x>=xco);
    eps_waypoints=xco;
    %Compute safety
    [Baseline_Grid(i,5),Baseline_Grid(i,6)]=wrap_Puma_RefSection(Ne,Nd,Baseline_Grid(i,1),f_eps,eps_waypoints,solparams);
end

toc
%% Plot results

if any(Baseline_Grid(i,6)==false)
    disp('The algorithm has not converged for some grid points')
else
    disp('Convergence everywhere')
end

figure(1)
clf
axis([24,40,0,100])
hold on
for i=1:(NKtheta*Ncg)
    if Baseline_Grid(i,5)
        scatter(100*Baseline_Grid(i,2)/c,100*Baseline_Grid(i,1)/Ktheta_nom,'g*')
    else
        scatter(100*Baseline_Grid(i,2)/c,100*Baseline_Grid(i,1)/Ktheta_nom,'r*')
    end
end
xlabel('Blade CoG [\% Chord]','Interpreter','latex')
ylabel('Control Chain Stiffness [\% Nominal]','Interpreter','latex')
