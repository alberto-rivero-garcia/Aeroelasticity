function UF=CompressiblePKGoland(Ne,Nd,Umax,Nu)
%This function returns a flutter velocity for the Goland Wing, using Ne
%elements, Nd modes and a velocity sweep up to Umax (in m/s) with Nu
%points. It's adapted from the PK_Goland.m script, rearranging functions
%and making the result easier to use. It employs a Prandtl-Glauert
%correction for Clalpha when applying theodorsen aerodynamics. Since it is
%a scalar correction, uniform over the whole wing, I preintegrate the
%incompressible matrices and then correct them inside of the loop. 
% Define relevant parameters: 
L=6.096;        %[m]
c=1.829;        %[m]

%For aerodynamic model: 
b=c/2;              %[m]
e=(0.33*c-b)/b;     %[-]
rho=1.225;          %[kg/m^3]
Clalpha_inc=2*pi;        %[1/rad] Incompressible value
c_sl=sqrt(1.4*101325/1.225); %[m/s] speed of sound at sea level on an ISA standard day. 

% First, solve the coupled bending-torsion problem to determine modal shapes. 
[fr,qr,FEMdata]=FEM_CoupledBendingTorsion(Ne,Nd);

%Rearrange modal results for easier manipulation later on: 
modes_cell=cell(1,Nd);
for i=1:Nd
    modes_cell{i}=@(y)[EvaluateCoupledBendingDisplacement(y,qr(:,i),FEMdata.deltay,FEMdata.Connect);...
    EvaluateCoupledTorsionDisplacement(y,qr(:,i),FEMdata.deltay,FEMdata.Connect)];
end

% Proceed with PK analysis: 
%Modal Mass and Stiffness matrices: based on mass normalized vectors. 
MM=eye(Nd,Nd);
KK=diag(fr.^2);

% Prepare aerodynamic matrices: preintegration is done to improve performance
%Define basic matrices
Ma=rho*b^2*[pi,-pi*b*e;pi*b*e,-pi*b^2*(1/8+e^2)];
Ca=rho*b^2*[0,pi;0,-pi*b*(1/2-e)];
Cw=[0,1];
Cwhat=[1,b*(1/2-e)];
Bw_inc=rho*b*Clalpha_inc*[1;b*(e+1/2)]; %incompressible value

%Convert into modal coordinates: 
Kstar1=nan(Nd,Nd);
Kstar2_inc=nan(Nd,Nd);
Kstar3_inc=nan(Nd,Nd);
Cstar1=nan(Nd,Nd);
Cstar2_inc=nan(Nd,Nd);
Cstar3_inc=nan(Nd,Nd);

%Use a for loop for integrals, leveraging the previously defined cell. 
for i=1:Nd
    for j=1:Nd
        %For aerodynamic stiffness matrix:
        Kstar1(i,j)=integral(@(y)AeroIntegrand(y,modes_cell{i},modes_cell{j},Ma),0,L);
        Kstar2_inc(i,j)=integral(@(y)AeroIntegrand(y,modes_cell{i},modes_cell{j},Bw_inc*Cw),0,L);
        Kstar3_inc(i,j)=integral(@(y)AeroIntegrand(y,modes_cell{i},modes_cell{j},Bw_inc*Cwhat),0,L);

        %For aerodynamic damping matrix:
        Cstar1(i,j)=integral(@(y)AeroIntegrand(y,modes_cell{i},modes_cell{j},Ca),0,L);
        Cstar2_inc(i,j)=integral(@(y)AeroIntegrand(y,modes_cell{i},modes_cell{j},Bw_inc*Cw),0,L);
        Cstar3_inc(i,j)=integral(@(y)AeroIntegrand(y,modes_cell{i},modes_cell{j},Bw_inc*Cwhat),0,L);

    end
end

% Set up pk loop: 

%Choose velocity, note we start from in vacuo response to use our knowledge
%of the structural system. 
UU=linspace(0,Umax,Nu); %[m/s]

%Define and Initialize eigenvalue and eigenvector storage matrices: 
Eigenvalues=nan(Nd,length(UU));
qq=nan(Nd,Nd,length(UU));

%Initialized with in vacuo response 
Eigenvalues(:,1)=fr*1i;
qq(:,:,1)=eye(Nd,Nd); %since we start with proper orthogonal modes

%Perform PK loop:
tol=1e-6;   %error tolerance in iterative solution
maxiter=20; %maximum number of iterations for each mode and velocity -> needs more than last time due to the more complicated equations (I'm using the secant method for the first initialization)
UF=nan; %easy to identify if no flutter is found.  
for i=2:length(UU)
    U=UU(i);
    Um1=UU(i-1);
    Beta=sqrt(1-(U/c_sl)^2); %for Prandtl-Glauert Correction
    Kstar2=1/Beta*Kstar2_inc;
    Kstar3=1/Beta*Kstar3_inc;
    Cstar2=1/Beta*Cstar2_inc;
    Cstar3=1/Beta*Cstar3_inc;
    for j=1:Nd
        %Prepare first guess for iteration: 
        lambda=Eigenvalues(j,i-1);
        q=qq(:,j,i-1);
        if Um1==0 %to catch issues with the first iteration, a bit of a hotfix, to be improved in the future
            Ctheo=zeros(Nd,Nd);
            Ktheo=zeros(Nd,Nd);
        else
            k=imag(lambda)*b/Um1;
            Ck=besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k));
            Ktheo=-(k/b)^2*Kstar1+real(Ck)*Kstar2-k/b*imag(Ck)*Kstar3;
            Ctheo=Cstar1+b/k*imag(Ck)*Cstar2+real(Ck)*Cstar3;
        end
        ddu=[2*MM*lambda*q-Um1*Ctheo*q, MM*lambda^2-Um1*Ctheo*lambda+(KK-Um1^2*Ktheo);0,2*q']\[Ctheo*lambda*q-(KK-2*Um1*Ktheo)*q;0];
        lambda=lambda+ddu(1)*(U-Um1);
        q=q+ddu(2:end)*(U-Um1);
        %Check after initial update: 
        k=imag(lambda)*b/U;
        Ck=besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k));
        Ktheo=-(k/b)^2*Kstar1+real(Ck)*Kstar2-k/b*imag(Ck)*Kstar3;
        Ctheo=Cstar1+b/k*imag(Ck)*Cstar2+real(Ck)*Cstar3;
        error=norm([(MM*lambda^2-U*Ctheo*lambda+(KK-U^2*Ktheo))*q;1-q'*q]);

        %Perform iteration until converged or maxed out iterations:
        iter=0;
        while error>tol && iter<=maxiter
            %Compute delta-lambda and delta-q
            A=[2*lambda*MM*q-U*Ctheo*q,MM*lambda^2-U*Ctheo*lambda+(KK-U^2*Ktheo);0,2*q'];
            v=[-MM*lambda^2*q+U*Ctheo*lambda*q-(KK-U^2*Ktheo)*q;1-q'*q];
            delta=A\v;
            %Update lambda and q
            lambda=lambda+delta(1);
            q=q+delta(2:end);
            %Update aerodynamic matrices:
            k=imag(lambda)*b/U;
            Ck=besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k));
            Ktheo=-(k/b)^2*Kstar1+real(Ck)*Kstar2-k/b*imag(Ck)*Kstar3;
            Ctheo=Cstar1+b/k*imag(Ck)*Cstar2+real(Ck)*Cstar3;
            %Update loop:
            iter=iter+1;
            %Recompute error: 
            error=norm([(MM*lambda^2-U*Ctheo*lambda+(KK-U^2*Ktheo))*q;1-q'*q]);
        end
        %Store new values:
        Eigenvalues(j,i)=lambda;
        qq(:,j,i)=q;
    end
    if max(real(Eigenvalues(:,i)))>0 %there's been flutter at some mode. 
        %UF=(U+Um1)/2; %simple interpolation, can be further refined if needed.
        %Advanced interpolation: 
        [~,I]=max(real(Eigenvalues(:,i)));
        gammaim1=real(Eigenvalues(I,i-1));
        gammai=real(Eigenvalues(I,i));
        UF=Um1-gammaim1/(gammai-gammaim1)*(U-Um1);
        break
    end
end

end

