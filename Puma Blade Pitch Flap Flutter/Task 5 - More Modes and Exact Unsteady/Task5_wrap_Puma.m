function [safe,convergence]=Task5_wrap_Puma(Ne,Nd,Ktheta,f_eps,eps_waypoints,solparams)
%This function wrapper is devised to handle the aeroelastic analysis of a
%given point in the stability map. It will simply return a value of true if
%the point is safe and a value of false if it is not safe. 

%In Task 5 it is modified to account for the exact distribution of reduced
%frequencies at hover. 

%It is based on v5 of the reference section wrapper, which uses the 
% continuation algorithm developed in v3, modified to handle
%cases where the in-vacuo system is already unstable. (statically
%unstable). These are automatically deemed unstable. It also uses v3 of the
%FEM, which ensures a PD mass matrix like v2(key!), but also checks the
%structure of the eigenvalues to catch possible in vacuo instabilities.
%Finally, it has reworked the treatment of cases with k close to 0, since
%the previous structure didn't work as intended. (It never triggered
%actually).

 

Nomega=solparams.Nomega;
maxiter=solparams.maxiter;
tol=solparams.tol;

%Physical parameters: 
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
T=[-1,0;0,1]; %for transformation of aerodynamic sub-matrices from h to w. 


%Perform modal analysis: 
[fr,qr,FEMdata,FEMCheck]=Task5_Parametric_FEM_RotatingCoupledBendingTorsion(Ne,Nd,Ktheta,f_eps,eps_waypoints);


if ~FEMCheck %if system doesn't meet sanity check in FEM
    safe=0;
    convergence=1;
else %if system is statically stable 
    %Rearrange modal results for easier manipulation later on: 
    modes_cell=cell(1,Nd);
    for i=1:Nd
        modes_cell{i}=@(y)[Offset_EvaluateCoupledBendingDisplacement(y,qr(:,i),FEMdata.deltax,FEMdata.Connect,FEMdata.xmin);...
        Offset_EvaluateCoupledTorsionDisplacement(y,qr(:,i),FEMdata.deltax,FEMdata.Connect,FEMdata.xmin)];
    end
    
    %Leverage modal solution for initial matrices: 
    MM=eye(Nd,Nd);
    KK=diag(fr.^2);
    
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
                %Construct k(x)
                k=@(x)imag(lambda_i)*b/Omega_v(i-1)./x;
                %Integrate Ktheohat and Ctheohat
                Ktheohat=nan(Nd,Nd);
                Ctheohat=nan(Nd,Nd);
                Kaer=@(x) -(k(x)/b)^2*Ma-k(x)/b*imag(Ck(k(x)))*Bw*Cwhat+real(Ck(k(x)))*Bw*Cw;
                Caer=@(x) Ca+real(Ck(k(x)))*Bw*Cwhat+b/k(x)*imag(Ck(k(x)))*Bw*Cw;
                for l=1:Nd
                    for m=1:Nd
                        Ktheohat(l,m)=integral(@(x)modes_cell{l}(x)'*x^2*Kaer(x)*T*modes_cell{m}(x),xco,R,'ArrayValued',true);
                        Ctheohat(l,m)=integral(@(x)modes_cell{l}(x)'*x*Caer(x)*T*modes_cell{m}(x),xco,R,'ArrayValued',true);
                    end
                end
            end 
            %Propagate: 
            DA=[2*MM*q_i*lambda_i-Omega_v(i-1)*Ctheohat*q_i,...
                MM*lambda_i^2-Omega_v(i-1)*Ctheohat*lambda_i+KK-Omega_v(i-1)^2*Ktheohat;0,2*q_i'];
            Db=[Ctheohat*lambda_i*q_i+2*Omega_v(i-1)*Ktheohat*q_i;0];
            ddOmega=DA\Db;
            lambda_i=lambda_i+ddOmega(1)*(Omega_v(i)-Omega_v(i-1));
            q_i=q_i+ddOmega(2:end)*(Omega_v(i)-Omega_v(i-1));
            %Update aerodynamic matrices: to new Omega and new prediction
            %Construct k(x)
            k=@(x)imag(lambda_i)*b/Omega_v(i)./x;
            %Integrate Ktheohat and Ctheohat
            Ktheohat=nan(Nd,Nd);
            Ctheohat=nan(Nd,Nd);
            Kaer=@(x) -(k(x)/b)^2*Ma-k(x)/b*imag(Ck(k(x)))*Bw*Cwhat+real(Ck(k(x)))*Bw*Cw;
            Caer=@(x) Ca+real(Ck(k(x)))*Bw*Cwhat+b/k(x)*imag(Ck(k(x)))*Bw*Cw;
            for l=1:Nd
                for m=1:Nd
                    Ktheohat(l,m)=integral(@(x)modes_cell{l}(x)'*x^2*Kaer(x)*T*modes_cell{m}(x),xco,R,'ArrayValued',true);
                    Ctheohat(l,m)=integral(@(x)modes_cell{l}(x)'*x*Caer(x)*T*modes_cell{m}(x),xco,R,'ArrayValued',true);
                end
            end
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
                %Construct k(x)
                k=@(x)imag(lambda_i)*b/Omega_v(i)./x;
                %Integrate Ktheohat and Ctheohat
                Ktheohat=nan(Nd,Nd);
                Ctheohat=nan(Nd,Nd);
                Kaer=@(x) -(k(x)/b)^2*Ma-k(x)/b*imag(Ck(k(x)))*Bw*Cwhat+real(Ck(k(x)))*Bw*Cw;
                Caer=@(x) Ca+real(Ck(k(x)))*Bw*Cwhat+b/k(x)*imag(Ck(k(x)))*Bw*Cw;
                for l=1:Nd
                    for m=1:Nd
                        Ktheohat(l,m)=integral(@(x)modes_cell{l}(x)'*x^2*Kaer(x)*T*modes_cell{m}(x),xco,R,'ArrayValued',true);
                        Ctheohat(l,m)=integral(@(x)modes_cell{l}(x)'*x*Caer(x)*T*modes_cell{m}(x),xco,R,'ArrayValued',true);
                    end
                end
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
    
    
    safe=~any(real(lambdas)>0); %unsafe if some eigenvalue has positive real part. 
    convergence=~any(conv==false); %converged if all conv flags are true. I only care about the last step for this.
end
end

function Ck=Ck(k)
%This function will evaluate Ck, with the limiting trick used to help
%convergence in divergence scenarios. It is designed to handle vector
%inputs just in case. 
Ck=nan(size(k));
for i=1:length(k)
    if k(i)<1e-6
        Ck(i)=1;
    else
        Ck(i)=besselh(1,2,k(i))/(besselh(1,2,k(i))+1i*besselh(0,2,k(i)));
    end
end

end