function UF=FlapPKGoland(Ne,Nd,Umax,Nu,Display)
%This function returns a flutter velocity for the Goland Wing, using Ne
%elements, Nd modes and a velocity sweep up to Umax (in m/s) with Nu
%points. It's adapted from the PK_Goland.m script, rearranging functions
%and making the result easier to use. Furthermore, it introduces the flap
%into the formulation to have a preliminary aeroservoelastic assessment.
%The Display argument allows the plotting of results. 
% Define relevant parameters: 

L=6.096;        %[m]
c=1.829;        %[m]
cf=0.45;        %[m]
%For aerodynamic model: 
b=c/2;              %[m]
ehat=(0.33*c-b)/b;  %[-]
chat=(c-cf-b)/b;    %[-]
rho=1.225;          %[kg/m^3]
Clalpha=2*pi;       %[1/rad]
yL=3.35;
yU=3.35+1.8;
%Define further aerodynamic models: 
mu=acos(chat);
T1=-1/3*sqrt(1-chat^2)*(2+chat^2)+chat*mu;
T2=chat*(1-chat^2)-sqrt(1-chat^2)*(1+chat^2)*mu+chat*mu^2;
T3=-(1/8+chat^2)*mu^2+1/4*chat*sqrt(1-chat^2)*mu*(7+2*chat^2)-1/8*(1-chat^2)*(5*chat^2+4);
T4=-mu+chat*sqrt(1-chat^2);
T5=-(1-chat^2)-mu^2+2*chat*sqrt(1-chat^2)*mu;
T7=-(1/8+chat^2)*mu+1/8*chat*sqrt(1-chat^2)*(7+2*chat^2);
T8=-1/3*sqrt(1-chat^2)*(2*chat^2+1)+chat*mu;
T9=1/2*(1/3*sqrt(1-chat^2)^3+ehat*T4);
T10=sqrt(1-chat^2)+mu;
T11=mu*(1-2*chat)+sqrt(1-chat^2)*(2-chat);
T12=sqrt(1-chat^2)*(2+chat)-mu*(2*chat+1);
T13=-1/2*(T7+(chat-ehat)*T1);


% First, solve the coupled bending-torsion problem to determine modal shapes. 
[fr,qr,FEMdata]=FEM_FlapCoupledBendingTorsion(Ne,Nd);

%Rearrange modal results for easier manipulation later on: 
modes_cell=cell(1,Nd);
for i=1:Nd
    modes_cell{i}=@(y)[EvaluateCoupledBendingDisplacement(y,qr(:,i),FEMdata.deltay,FEMdata.Connect);...
    EvaluateCoupledTorsionDisplacement(y,qr(:,i),FEMdata.deltay,FEMdata.Connect);Betafun(y,qr(:,i))];
end

% Proceed with PK analysis: 
%Modal Mass and Stiffness matrices: based on mass normalized vectors. 
MM=eye(Nd,Nd);
KK=diag(fr.^2);

% Prepare aerodynamic matrices: preintegration is done to improve performance
%Define basic matrices
Ca=rho*b^2*[0,pi,-T4;0,-pi*b*(1/2-ehat),-b*(T1-T8-(chat-ehat)*T4+T11/2);...
    0,-b*(-2*T9-T1+T4*(ehat-1/2)),b*T4*T11/2/pi];
Ka=rho*b^2*[0,0,0;0,0,-(T4+T10);0,0,-(T5-T4*T10)/pi];
Ma=rho*b^2*[pi,-pi*b*ehat,-b*T1;pi*b*ehat,-pi*b^2*(1/8+ehat^2),b^2*(T7+(chat-ehat)*T1);...
    b*T1,-2*b^2*T13,b^2*T3/pi];
Bw=rho*b*Clalpha*[1;b*(ehat+1/2);-T12*b/2/pi];
Cw=[0,1,T10/pi];
Cwhat=[1,b*(1/2-ehat),b*T11/2/pi];

%Convert into modal coordinates: 
Kstar1=nan(Nd,Nd);
Kstar2=nan(Nd,Nd);
Kstar3=nan(Nd,Nd);
Kstar4=nan(Nd,Nd);
Cstar1=nan(Nd,Nd);
Cstar2=nan(Nd,Nd);
Cstar3=nan(Nd,Nd);

%Use a for loop for integrals, leveraging the previously defined cell. Needs some rework to adapt to the new Kaer and Caer models.  
for i=1:Nd
    for j=1:Nd
        %For aerodynamic stiffness matrix:
        Kstar1(i,j)=integral(@(y)AeroIntegrand(y,modes_cell{i},modes_cell{j},Ma),0,L,'Waypoints',[yL,yU]);
        Kstar2(i,j)=integral(@(y)AeroIntegrand(y,modes_cell{i},modes_cell{j},Bw*Cw),0,L,'Waypoints',[yL,yU]);
        Kstar3(i,j)=integral(@(y)AeroIntegrand(y,modes_cell{i},modes_cell{j},Bw*Cwhat),0,L,'Waypoints',[yL,yU]);
        Kstar4(i,j)=integral(@(y)AeroIntegrand(y,modes_cell{i},modes_cell{j},Ka),0,L,'Waypoints',[yL,yU]);

        %For aerodynamic damping matrix:
        Cstar1(i,j)=integral(@(y)AeroIntegrand(y,modes_cell{i},modes_cell{j},Ca),0,L,'Waypoints',[yL,yU]);
        Cstar2(i,j)=integral(@(y)AeroIntegrand(y,modes_cell{i},modes_cell{j},Bw*Cw),0,L,'Waypoints',[yL,yU]);
        Cstar3(i,j)=integral(@(y)AeroIntegrand(y,modes_cell{i},modes_cell{j},Bw*Cwhat),0,L,'Waypoints',[yL,yU]);

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
            Ktheo=-(k/b)^2*Kstar1+real(Ck)*Kstar2-k/b*imag(Ck)*Kstar3+Kstar4;
            Ctheo=Cstar1+b/k*imag(Ck)*Cstar2+real(Ck)*Cstar3;
        end
        ddu=[2*MM*lambda*q-Um1*Ctheo*q, MM*lambda^2-Um1*Ctheo*lambda+(KK-Um1^2*Ktheo);0,2*q']\[Ctheo*lambda*q-(KK-2*Um1*Ktheo)*q;0];
        lambda=lambda+ddu(1)*(U-Um1);
        q=q+ddu(2:end)*(U-Um1);
        %Check after initial update: 
        k=imag(lambda)*b/U;
        Ck=besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k));
        Ktheo=-(k/b)^2*Kstar1+real(Ck)*Kstar2-k/b*imag(Ck)*Kstar3+Kstar4;
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
            Ktheo=-(k/b)^2*Kstar1+real(Ck)*Kstar2-k/b*imag(Ck)*Kstar3+Kstar4;
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
        iF=i;
        break
    end
end

if Display
    %Plot frequencies
    Frequencies=imag(Eigenvalues)/2/pi; %[Hz]
    
    figure 
    for i=1:Nd
        curvelabel='Mode' + string(i);
        plot(UU(1:iF),Frequencies(i,1:iF),'DisplayName',curvelabel)
        hold on
    end
    grid minor
    legend
    title('Modal Frequencies vs Airspeed')
    xlabel('$U \; [m/s]$','Interpreter','latex')
    ylabel('$f \; [Hz]$','Interpreter','latex')

    %Another option is simply the real part of the eigenvalue: 
    figure
    for i=1:Nd
        curvelabel="Mode"+i;
        plot(UU(1:iF),real(Eigenvalues(i,1:iF)),'DisplayName',curvelabel)
        hold on
    end
    grid minor 
    legend
    title('Modal real part vs Airspeed')
    xlabel('$U \; [m/s]$','Interpreter','latex')
    ylabel('$real(\lambda) \; [1/s]$','Interpreter','latex')

end

end

