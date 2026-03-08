function [KK_local,MM_local]=Parametric_RotatingLocalElementMatrices(xl,xr,Omega,R,f_eps,eps_waypoints)
%This function computes the local stiffness and mass matrices of the
%element contained between xl and xr of the Puma blade. It relies on
%numerical integration of my development of the PVW and my version of the
%input data. The blades are assumed to be rotating at Omega rad/s.

%This function has been modified to accept a displacement epsilon of the
%section cg's, in order to compute the stability maps. The displacement is
%given as a function f_eps and eps_waypoints are included in the relevant
%integration. 


%Preallocate matrices:
KK_local=zeros(7,7);
MM_local=zeros(7,7);

%Define xi and all of the corresponding shape functions and their relevant
%derivatives. Arrange them in cells afterwards to really ease programming: 
xi=@(x)(x-xl)/(xr-xl);

psiw1=@(x) 1-3*xi(x).^2+2*xi(x).^3;
psiw2=@(x)(xi(x)-2*xi(x).^2+xi(x).^3)*(xr-xl);
psiw3=@(x) 3*xi(x).^2-2*xi(x).^3;
psiw4=@(x) (-xi(x).^2+xi(x).^3)*(xr-xl);
psit1=@(x) 1-3*xi(x)+2*xi(x).^2;
psit2=@(x) 4*xi(x)-4*xi(x).^2;
psit3=@(x) 2*xi(x).^2-xi(x);

psi_cell={psiw1,psiw2,psiw3,psiw4,psit1,psit2,psit3};

dpsiw1dx=@(x)(-6*xi(x)+6*xi(x).^2)/(xr-xl);
dpsiw2dx=@(x)(1-4*xi(x)+3*xi(x).^2);
dpsiw3dx=@(x)(6*xi(x)-6*xi(x).^2)/(xr-xl);
dpsiw4dx=@(x)(-2*xi(x)+3*xi(x).^2);
dpsit1dx=@(x)(-3+4*xi(x))/(xr-xl);
dpsit2dx=@(x)(4-8*xi(x))/(xr-xl);
dpsit3dx=@(x)(4*xi(x)-1)/(xr-xl);

dpsidx_cell={dpsiw1dx,dpsiw2dx,dpsiw3dx,dpsiw4dx,dpsit1dx,dpsit2dx,dpsit3dx};

d2psiw1dx2=@(x)(-6+12*xi(x))/((xr-xl)^2);
d2psiw2dx2=@(x)(-4+6*xi(x))/(xr-xl);
d2psiw3dx2=@(x)(6-12*xi(x))/((xr-xl)^2);
d2psiw4dx2=@(x)(-2+6*xi(x))/(xr-xl);

d2psidx2_cell={d2psiw1dx2,d2psiw2dx2,d2psiw3dx2,d2psiw4dx2};

%Integrate term by term of the PVW accordingly and add them to the matrices:
%The first call to the data functions identifies waypoints that can be
%provided to the integral function. Note that if none are found, it will
%return an empty vector, which is simply ignored by the integral function.


%Stiffness matrix: 
%Classic bending stiffness: 
[~,waypoints]=EI_Puma([xl,xr]);
for i=1:4
    for j=1:4
        KK_local(i,j)=KK_local(i,j)+integral(@(x)d2psidx2_cell{i}(x).*EI_Puma(x).*d2psidx2_cell{j}(x),xl,xr,"Waypoints",waypoints);
    end
end

%Centrifugal bending stiffness: 

%Construct an in-line function that will compute N(x)
%Proceed with element integration: N(x) will not be differentiable in them,
%so I use them just in case. N(x) defined as another function to be able to
%handle vector inputs in x. 

[~,waypoints]=m_Puma([xl,xr]);
for i=1:4
    for j=1:4
        KK_local(i,j)=KK_local(i,j)+integral(@(x)dpsidx_cell{i}(x).*N(x,Omega,R).*dpsidx_cell{j}(x),xl,xr,'Waypoints',waypoints);
    end
end

%Classic torsional stiffness: 
[~,waypoints]=GJ_Puma([xl,xr]);
for i=5:7
    for j=5:7
        KK_local(i,j)=KK_local(i,j)+integral(@(x)dpsidx_cell{i}(x).*GJ_Puma(x).*dpsidx_cell{j}(x),xl,xr,"Waypoints",waypoints);
    end
end

%Torsion effect on bending dofs: requires combining waypoints
[~,waypoints_1]=m_Puma([xl,xr]);
[~,waypoints_2]=ycg_Puma([xl,xr]);
waypoints=sort(unique([waypoints_1;waypoints_2;eps_waypoints])); %coding implementation of the logic checked. 
check=zeros(7,7);
for i=1:4
    for j=5:7
        KK_local(i,j)=KK_local(i,j)+Omega^2*integral(@(x)dpsidx_cell{i}(x).*m_Puma(x).*(ycg_Puma(x)+f_eps(x)).*x.*psi_cell{j}(x),xl,xr,'Waypoints',waypoints);
        check(i,j)=check(i,j)+Omega^2*integral(@(x)dpsidx_cell{i}(x).*m_Puma(x).*(ycg_Puma(x)+f_eps(x)).*x.*psi_cell{j}(x),xl,xr,'Waypoints',waypoints);
    end
end

%Bending effect on torsion dofs: can reuse the above waypoints, but I'll
%rewrite it for robustness (it's a 1 time code, no need for performance)
[~,waypoints_1]=m_Puma([xl,xr]);
[~,waypoints_2]=ycg_Puma([xl,xr]);
waypoints=sort(unique([waypoints_1;waypoints_2;eps_waypoints])); %coding implementation of the logic checked. 

for i=5:7
    for j=1:4
        KK_local(i,j)=KK_local(i,j)+Omega^2*integral(@(x)psi_cell{i}(x).*m_Puma(x).*(ycg_Puma(x)+f_eps(x)).*x.*dpsidx_cell{j}(x),xl,xr,'Waypoints',waypoints);
        check(i,j)=check(i,j)+Omega^2*integral(@(x)psi_cell{i}(x).*m_Puma(x).*(ycg_Puma(x)+f_eps(x)).*x.*dpsidx_cell{j}(x),xl,xr,'Waypoints',waypoints);
    end
end

%Propeller moment: 
[~,waypoints_1]=Itheta_Puma([xl,xr]);
[~,waypoints_2]=m_Puma([xl,xr]);
[~,waypoints_3]=ycg_Puma([xl,xr]);
waypoints=sort(unique([waypoints_1;waypoints_2;waypoints_3;eps_waypoints]));

for i=5:7
    for j=5:7
        KK_local(i,j)=KK_local(i,j)+Omega^2*integral(@(x)psi_cell{i}(x).*(Itheta_Puma(x)+m_Puma(x).*(f_eps(x).^2+2*f_eps(x).*ycg_Puma(x))).*psi_cell{j}(x),xl,xr,'Waypoints',waypoints);
    end
end

%Mass matrix: 
%Classic w on w: 
[~,waypoints]=m_Puma([xl,xr]);
for i=1:4
    for j=1:4
        MM_local(i,j)=MM_local(i,j)+integral(@(x)psi_cell{i}(x).*m_Puma(x).*psi_cell{j}(x),xl,xr,'Waypoints',waypoints);
    end
end

%Classic theta on w: 
[~,waypoints_1]=m_Puma([xl,xr]);
[~,waypoints_2]=ycg_Puma([xl,xr]);
waypoints=sort(unique([waypoints_1;waypoints_2;eps_waypoints])); %coding implementation of the logic checked. 

for i=1:4
    for j=5:7
        MM_local(i,j)=MM_local(i,j)+integral(@(x)psi_cell{i}(x).*m_Puma(x).*(ycg_Puma(x)+f_eps(x)).*psi_cell{j}(x),xl,xr,'Waypoints',waypoints);
    end
end

%Classic w on theta:
[~,waypoints_1]=m_Puma([xl,xr]);
[~,waypoints_2]=ycg_Puma([xl,xr]);
waypoints=sort(unique([waypoints_1;waypoints_2;eps_waypoints])); %coding implementation of the logic checked. 

for i=5:7
    for j=1:4
        MM_local(i,j)=MM_local(i,j)+integral(@(x)psi_cell{i}(x).*m_Puma(x).*(ycg_Puma(x)+f_eps(x)).*psi_cell{j}(x),xl,xr,'Waypoints',waypoints);
    end
end

%Classic theta on theta: 
[~,waypoints_1]=Itheta_Puma([xl,xr]);
[~,waypoints_2]=m_Puma([xl,xr]);
[~,waypoints_3]=ycg_Puma([xl,xr]);
waypoints=sort(unique([waypoints_1;waypoints_2;waypoints_3;eps_waypoints]));

for i=5:7
    for j=5:7
        MM_local(i,j)=MM_local(i,j)+integral(@(x)psi_cell{i}(x).*(Itheta_Puma(x)+m_Puma(x).*(f_eps(x).^2+2*f_eps(x).*ycg_Puma(x))).*psi_cell{j}(x),xl,xr,'Waypoints',waypoints);
    end
end
end

function N_x=N(x,Omega,R)
N_x=nan(size(x));
for i=1:length(x)
    [~,waypoints]=m_Puma([x(i),R]);
    N_x(i)=Omega^2*integral(@(xprime)m_Puma(xprime).*xprime,x(i),R,'Waypoints',waypoints);
end
end