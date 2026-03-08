%This script analyzes the aeroelastic stability of the Puma's blade in
%hover, using a pseudo-Theodorsen model for aerodynamics. It employs
%Theodorsen's expressions and strip theory, except for C(k), which is taken
%at a reference section at 0.75R.


%Basic file management
addpath("Data\","FEM\")
tic
%Define solution parameters: 
Ne=20;
Nd=2;  %tbd after we choose the first torsion and bending modes
Nomega=20;
maxiter=20;
tol=1e-8;

%Define problem parameters: 
xe= 0.289;      %[m]
xco=2.080;      %[m]
R=7.490;        %[m]
c=0.537;        %[m]
Omega=4.5*2*pi; %[rad/s]

%For aerodynamic model:
b=c/2;          %[m]
Clalpha=2*pi;   %[-]
rho=1.225;      %[kg/m^3]

%Construct baseline aerodynamic matrices: 
Ma=pi*rho*b^2*[1,b/2;-b/2,-3/8*b^2];
Ca=pi*rho*b^2*[0,1;0,-b];
Bw=rho*b*Clalpha*[1;0];
Cw=[0,1];
Cwhat=[1,b];


%Extract desired modes using FEM:
[fr,qr,FEMdata]=FEM_RotatingCoupledBendingTorsion(Ne,Nd);

%Rearrange modal results for easier manipulation later on: 
modes_cell=cell(1,Nd);
for i=1:Nd
    modes_cell{i}=@(y)[Offset_EvaluateCoupledBendingDisplacement(y,qr(:,i),FEMdata.deltax,FEMdata.Connect,FEMdata.xmin);...
    Offset_EvaluateCoupledTorsionDisplacement(y,qr(:,i),FEMdata.deltax,FEMdata.Connect,FEMdata.xmin)];
end

%Leverage modal solution for initial matrices: 
MM=eye(Nd,Nd);
KK=diag(fr.^2);

%Pre-integrate aerodynamic sub-matrices: 
Kstar1=nan(Nd,Nd);
Kstar2=nan(Nd,Nd);
Kstar3=nan(Nd,Nd);
Cstar1=nan(Nd,Nd);
Cstar2=nan(Nd,Nd);
Cstar3=nan(Nd,Nd);
T=[-1,0;0,1]; %for transformation of aerodynamic sub-matrices from h to w. 
for i=1:Nd
    for j=1:Nd
        Kstar1(i,j)=integral(@(x)modes_cell{i}(x)'*x^2*Ma*T*modes_cell{j}(x),xco,R,"ArrayValued",true);
        Kstar2(i,j)=integral(@(x)modes_cell{i}(x)'*x^2*Bw*Cwhat*T*modes_cell{j}(x),xco,R,"ArrayValued",true);
        Kstar3(i,j)=integral(@(x)modes_cell{i}(x)'*x^2*Bw*Cw*T*modes_cell{j}(x),xco,R,'ArrayValued',true);
        Cstar1(i,j)=integral(@(x)modes_cell{i}(x)'*x*Ca*T*modes_cell{j}(x),xco,R,'ArrayValued',true);
        Cstar2(i,j)=integral(@(x)modes_cell{i}(x)'*x*Bw*Cwhat*T*modes_cell{j}(x),xco,R,'ArrayValued',true);
        Cstar3(i,j)=integral(@(x)modes_cell{i}(x)'*x*Bw*Cw*T*modes_cell{j}(x),xco,R,'ArrayValued',true);
    end
end

%Proceed with the continuation algorithm:
%Construct the Omega vector: 
Omega_v=linspace(0,Omega,Nomega);

%Eigenvalues and eigenvectors, at the end of the loop they will reach the
%desired value. 
%Prepare initial guess: elastomechanical system
lambdas=1i*(fr);
qs=eye(Nd,Nd);
conv=nan(Nd,1);
%Continuation loop itself: 
for i=2:length(Omega_v)
    for j=1:Nd
        %Search for j-th modal solution.
        %Initialize search: start from previous value
        lambda_i=lambdas(j);
        q_i=qs(:,j);
        %Construct aero matrices for the corresponding mode at previous Omega: 
        if i==2 %to catch the fact that Omega=0 leads to undefined k 
            Ktheohat=zeros(Nd,Nd);
            Ctheohat=zeros(Nd,Nd);
        else
            k=imag(lambda_i)*b/(3/4*Omega_v(i-1)*R);
            Ck=besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k));
            RCk=real(Ck);
            ICk=imag(Ck);
            Ktheohat=-(k/b)^2*Kstar1-(k/b)*ICk*Kstar2+RCk*Kstar3;
            Ctheohat=Cstar1+RCk*Cstar2+b/k*ICk*Cstar3; 
        end
        %Propagate: 
        DA=[2*MM*q_i*lambda_i-Omega_v(i-1)*Ctheohat*q_i,...
            MM*lambda_i^2-Omega_v(i-1)*Ctheohat*lambda_i+KK-Omega_v(i-1)^2*Ktheohat;0,2*q_i'];
        Db=[Ctheohat*lambda_i*q_i+2*Omega_v(i-1)*Ktheohat*q_i;0];
        ddOmega=DA\Db;
        lambda_i=lambda_i+ddOmega(1)*(Omega_v(i)-Omega_v(i-1));
        q_i=q_i+ddOmega(2:end)*(Omega_v(i)-Omega_v(i-1));
        %Update aerodynamic matrices: 
        k=imag(lambda_i)*b/(3/4*Omega_v(i)*R);
        Ck=besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k));
        RCk=real(Ck);
        ICk=imag(Ck);
        Ktheohat=-(k/b)^2*Kstar1-(k/b)*ICk*Kstar2+RCk*Kstar3;
        Ctheohat=Cstar1+RCk*Cstar2+b/k*ICk*Cstar3;        
        %Check initial error:
        e_v=[(MM*lambda_i^2-Omega_v(i)*Ctheohat*lambda_i+KK-Omega_v(i)^2*Ktheohat)*q_i;q_i'*q_i-1];
        %Perform while loop:
        iter=0;
        while norm(e_v)>tol && iter<=maxiter
            %Compute delta corrections: 
            A=[2*MM*lambda_i*q_i-Omega_v(i)*Ctheohat*q_i,...
                MM*lambda_i^2-Omega_v(i)*Ctheohat*lambda_i+KK-Omega_v(i)^2*Ktheohat;0,2*q_i'];
            delta_v=A\(-e_v);
            %Update guess: 
            lambda_i=lambda_i+delta_v(1);
            q_i=q_i+delta_v(2:end);
            %Recompute aerodynamic matrices: 
            k=imag(lambda_i)*b/(3/4*Omega_v(i)*R);
            Ck=besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k));
            RCk=real(Ck);
            ICk=imag(Ck);
            Ktheohat=-(k/b)^2*Kstar1-(k/b)*ICk*Kstar2+RCk*Kstar3;
            Ctheohat=Cstar1+RCk*Cstar2+b/k*ICk*Cstar3;
            %Recompute error: 
            e_v=[(MM*lambda_i^2-Omega_v(i)*Ctheohat*lambda_i+KK-Omega_v(i)^2*Ktheohat)*q_i;q_i'*q_i-1];
            %Update iteration counter: 
            iter=iter+1;
        end
        %Update results: 
        lambdas(j)=lambda_i;
        qs(:,j)=q_i;
        %Store convergence results: 
        if norm(e_v)<=tol
            conv(j)=true;
        else
            conv(j)=false;
        end
    end
end

if any(conv==false)
    disp('The method has not converged properly')
end

%Final postprocessing: 
display(lambdas)

damping_ratios=nan(size(lambdas));

for i=1:Nd
    damping_ratios(i)=-real(lambdas(i))/abs(lambdas(i));
end

display(damping_ratios)

toc