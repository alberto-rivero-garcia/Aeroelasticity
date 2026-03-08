%This script analyzes the aeroelastic stability of the Puma's blade in
%hover, using a Quasisteady model for aerodynamics. 

%Basic file management
addpath("Data\","FEM\")
tic
%Define problem parameters: 
Ne=20;
Nd=2;  

xe= 0.289;      %[m]
xco=2.080;      %[m]
R=7.490;        %[m]
c=0.537;        %[m]
Omega=4.5*2*pi; %[rad/s]

%For aerodynamic model:
b=c/2;          %[m]
Clalpha=2*pi;   %[-]
rho=1.225;      %[kg/m^3]
Ka=[0,rho*b*Clalpha;0,0];
Ca=rho*b*[-Clalpha,pi*b+b*Clalpha;0,-pi*b^2];
Ma=pi*rho*b^2*[-1,b/2;b/2,-3/8*b^2];


%Extract desired modes using FEM:
[fr,qr,FEMdata]=FEM_RotatingCoupledBendingTorsion(Ne,Nd);

%Rearrange modal results for easier manipulation later on: 
modes_cell=cell(1,Nd);
for i=1:Nd
    modes_cell{i}=@(y)[Offset_EvaluateCoupledBendingDisplacement(y,qr(:,i),FEMdata.deltax,FEMdata.Connect,FEMdata.xmin);...
    Offset_EvaluateCoupledTorsionDisplacement(y,qr(:,i),FEMdata.deltax,FEMdata.Connect,FEMdata.xmin)];
end

% Proceed with PK analysis: 
%Modal Mass and Stiffness matrices: based on mass normalized vectors. 
MM=eye(Nd,Nd);
KK=diag(fr.^2);

%Construct modal aerodynamic matrices: 
Kaer=nan(Nd,Nd);
Caer=nan(Nd,Nd);
Maer=nan(Nd,Nd);

for i=1:Nd
    for j=1:Nd       
        Kaer(i,j)=integral(@(x)(Omega*x)^2*modes_cell{i}(x)'*Ka*modes_cell{j}(x),xco,R,'ArrayValued',true);
        Caer(i,j)=integral(@(x)(Omega*x)*modes_cell{i}(x)'*Ca*modes_cell{j}(x),xco,R,'ArrayValued',true);
        Maer(i,j)=integral(@(x)modes_cell{i}(x)'*Ma*modes_cell{j}(x),xco,R,'ArrayValued',true);
    end
end

%Extract eigenvalues using Polyeig and display them: 
lambdas=polyeig(KK-Kaer,-Caer,MM-Maer);

display(lambdas)

%Extract roots with positive imaginary part and compute damping ratios: 
pos_lambdas=lambdas(imag(lambdas)>0);

damping_ratios=nan(size(pos_lambdas));

for i=1:length(pos_lambdas)
    damping_ratios(i)=-real(pos_lambdas(i))/abs(pos_lambdas(i));
end

display(damping_ratios)
toc