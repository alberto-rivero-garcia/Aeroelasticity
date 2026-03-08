%This script starts from a pre-computed baseline grid and refines it to get
%a better shape for the curve separating the safe and unsafe regions. It
%does so by performing a horizontal search method (so varying ycg at
%constant Ktheta). This choice is the result of the expected shape given by
%the class material, which would be messy to detect with a combined
%vertical shape. In any case, vertical refinement can be achieved by
%refining the original grid.

%The idea for the algorithm is: start at a line -> identify the change from
%safe to unsafe -> bisect ycg between these points -> perform pk analysis
%(based on reference section) -> add to grid.

%Main challenge: doing this repeatedly and programmatically -> so one can
%choose the desired level of refinement.

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%File management
addpath("Data\","FEM\")

%Set refinement level:
Nref=3; %set refinement level

%Specify solver parameters: 
Ne=20;
Nd=2;
solparams.Nomega=10;
solparams.maxiter=50;
solparams.tol=1e-6;

%Define physical parameters and prepare computation of epsilon:
c=0.537;            %[m]
Ktheta_nom=33032;   %[Nm/rad]
xe= 0.289;          %[m]
xco=2.080;          %[m]
R=7.490;            %[m]
[~,waypoints_1]=m_Puma([xe,R]);
[~,waypoints_2]=ycg_Puma([xe,R]);
waypoints=sort(unique([waypoints_1;waypoints_2]));
Blade_mass=integral(@(x)m_Puma(x),xe,R,'Waypoints',waypoints_1);
Blade_co_mass=integral(@(x)m_Puma(x),xco,R,'Waypoints',waypoints_1);
Ycg_baseline=integral(@(x)m_Puma(x).*ycg_Puma(x),xe,R,'Waypoints',waypoints)/Blade_mass;




%First step: load and search the baseline grid to add the corresponding
%points. 

load("Baseline_Grid.mat")

%Some basic grid sizing: 
B_NKtheta=length(unique(Baseline_Grid(:,1)));
B_NYcg=length(unique(Baseline_Grid(:,2)));
B_Npoints=B_NKtheta*B_NYcg;


%First refinement level: must always start with a full search.
%Search each row of constant NKtheta to find the change from safe to
%unsafe: 
Refined_Grid=[Baseline_Grid;nan(B_NKtheta,6)];
for i=1:B_NKtheta
    %Find change
    for j=1:(B_NYcg-1) %start from the left
        point=B_NYcg*(i-1)+j;
        right_point=Baseline_Grid(point,4);
        if ~Baseline_Grid(right_point,5) %if the point to the right is unsafe
            %Add a new point to the grid: midpoint between the others
            new_point_v=[Baseline_Grid(point,1),(Baseline_Grid(point,2)+Baseline_Grid(right_point,2))/2,...
                point,right_point,nan,nan];
            Refined_Grid(B_Npoints+i,:)=new_point_v;
            %Check the safety of the new point
            %Compute epsilon:
            Ycg_prime=c/4-Refined_Grid(B_Npoints+i,2);
            epsilon=(Ycg_prime-Ycg_baseline)*Blade_mass/Blade_co_mass;
            %Prepare cg input
            f_eps=@(x)epsilon*(x>=xco);
            eps_waypoints=xco;
            %Compute safety
            [Refined_Grid(B_Npoints+i,5),Refined_Grid(B_Npoints+i,6)]=v5_wrap_Puma_RefSection(Ne,Nd,Refined_Grid(B_Npoints+i,1),f_eps,eps_waypoints,solparams);
            %Break nested for loop (stop search in this row)
            break
        end
    end
end

%Remaining refinement levels: (these are quicker, since they only have to
%check the points added in the previous one: 
Refined_Grid=[Refined_Grid;nan((Nref-1)*B_NKtheta,6)];
for i=2:Nref
    %Go through each of the 10 added points: 
    for j=1:B_NKtheta
        point=B_Npoints+(i-2)*B_NKtheta+j;
        if Refined_Grid(point,5) %point is safe -> go to the right
        right_point=Refined_Grid(point,4);
        new_point_v=[Refined_Grid(point,1),(Refined_Grid(point,2)+Refined_Grid(right_point,2))/2,...
            point,right_point,nan,nan];
        else %point is unsafe -> go to the left
        left_point=Refined_Grid(point,3);
        new_point_v=[Refined_Grid(point,1),(Refined_Grid(point,2)+Refined_Grid(left_point,2))/2,...
            left_point,point,nan,nan];
        end
        new_point=point+B_NKtheta;
        Refined_Grid(new_point,:)=new_point_v;
        %Check the safety of the new point
        %Compute epsilon:
        Ycg_prime=c/4-Refined_Grid(new_point,2);
        epsilon=(Ycg_prime-Ycg_baseline)*Blade_mass/Blade_co_mass;
        %Prepare cg input
        f_eps=@(x)epsilon*(x>=xco);
        eps_waypoints=xco;
        %Compute safety
        [Refined_Grid(new_point,5),Refined_Grid(new_point,6)]=v5_wrap_Puma_RefSection(Ne,Nd,Refined_Grid(new_point,1),f_eps,eps_waypoints,solparams);
    end
end

%% 

if any(Refined_Grid(i,6)==false)
    disp('The algorithm has not converged for some grid points')
else
    disp('Convergence everywhere')
end

%Plot resulting grid: 
[n,~]=size(Refined_Grid);

figure
clf
axis([24,40,0,100])
hold on
for k=1:n
    if Refined_Grid(k,5)
        scatter(100*Refined_Grid(k,2)/c,100*Refined_Grid(k,1)/Ktheta_nom,'g*')
    else
        scatter(100*Refined_Grid(k,2)/c,100*Refined_Grid(k,1)/Ktheta_nom,'r*')
    end
end
xlabel('Blade CoG [\% Chord]','Interpreter','latex')
ylabel('Control Chain Stiffness [\% Nominal]','Interpreter','latex')