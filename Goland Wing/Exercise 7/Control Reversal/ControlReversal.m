%This script computes the control reversal speed for the Goland wing, using
%only torsion degrees of freedom. The procedure it employs is detailed in 
%the report. 
addpath("Auxiliary Functions\")
%Define number of elements for finite element discretization: 
Ne=15; 

%Physical Parameters: 
%Wing
L=6.096;        %[m]
c=1.829;        %[m]
shat=0.1;          %[-] Distance between elastic and mass axes, referred to the chord. 
EI=9.77*1e6;    %[N*m2]
m=35.72;        %[kg/m]
GJ=987600;      %[N m^2]
Io=7.452;       %[kg m] 
s=shat*c;
Itheta=Io+m*s^2;
xe=0.33*c;
%Flap:
Mf=8.92;            %[kg]
If=0.09;            %[kg m^2]
cf=0.45;            %[m]
yL=3.35;            %[m]
yU=yL+1.8;          %[m]
yf=(yL+yU)/2;       %[m]
xf=c-cf/2;          %[m]
Ih=If+(cf/2)^2*Mf;  %[kg m^2]
KBeta=6.48*1e3;     %[Nm/rad]


%Number of torsion nodes and dofs:
NnodesT=2*Ne+1;
NdofsT=NnodesT;
%Overall dofs: 
Ndofs=NdofsT+2; %now we must also consider Beta and Beta_0 as a dofs.

%Element size: 
deltay=L/Ne;

%Define connectivity matrix: 
 
Connect=nan(Ne,3);

for i=1:Ne
    Connect(i,:)=[1+2*(i-1),2+2*(i-1),3+2*(i-1)];
end


%Torsion problem:
K_local=GJ/deltay/3*[7,-8,1;-8,16,-8;1,-8,7];


%Assemble structural stiffness matrix: 
KK=zeros(Ndofs,Ndofs);
for i=1:Ne
    KK(Connect(i,:),Connect(i,:))=KK(Connect(i,:),Connect(i,:))+K_local;
end
KK(end-1,end-1)=KBeta;
KK(end-1,end)=-KBeta; %Beta_0 term

%Construct aerodynamic stiffness problem: 
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

%Define basic steady aero matrix: 
Kast=2*b*[0,Clalpha,Clalpha*T10/pi;0,b*(ehat+1/2)*Clalpha,b*((ehat+1/2)*Clalpha*T10/pi-(T4+T10));0,-b*T12/2/pi*Clalpha,-b*(T10*T12/2/pi^2*Clalpha+(T5-T4*T10)/pi)];

%Assemble global aerodynamic stiffness matrix. Each element requires
%integration of its local aero stiffness matrix. 
KKaer=zeros(Ndofs,Ndofs);
ymin=(0:(Ne-1))*deltay;
ymax=(1:Ne)*deltay;
for i=1:Ne
    Ka_local=TorsionLocalAeroMatrix(ymin(i),ymax(i),Kast(2:end,2:end));
    %For w and theta degrees of freedom:
    KKaer(Connect(i,:),Connect(i,:))=KKaer(Connect(i,:),Connect(i,:))+Ka_local(1:3,1:3);
    %For w and theta dofs, effect of beta:
    KKaer(Connect(i,:),end-1)=KKaer(Connect(i,:),end-1)+Ka_local(1:3,end);
    %For beta, effect of w and theta dofs: 
    KKaer(end-1,Connect(i,:))=KKaer(end-1,Connect(i,:))+Ka_local(end,1:3);
    %For beta:
    KKaer(end-1,end-1)=KKaer(end-1,end-1)+Ka_local(end,end);
end

%Add the lift equation: 
for i=1:Ne
    %Compute element lift contribution
    L_local=TorsionLiftContribution(ymin(i),ymax(i),Kast);
    %Add w and theta contribution:
    KKaer(end,Connect(i,:))=KKaer(end,Connect(i,:))+L_local(1:3);
    %Add beta contribution: 
    KKaer(end,end-1)=KKaer(end,end-1)+L_local(end);
end

%Apply degrees of freedom: take out those related to restricted dofs in
%both matrices: 
%Eliminate torsion of the root: 
KKr=KK(2:end,2:end);
KKaer_r=KKaer(2:end,2:end);

%Solve the eigenvalue problem: 
[Vr,D]=eig(KKr,KKaer_r);

[qv,sorting]=sort(real(diag(D)));
Vr=Vr(:,sorting);

%Extract information about control reversal, we want the lowest real, positive q with non-null beta term in the eigenvector. 
j=0;
qfound=0; %flag variable
qreversal=nan;
while j<length(qv) && qfound==0
    j=j+1;
    if imag(qv(j))==0
        if qv(j)>=1e-4 %tolerance to avoid the first analytical solution, which is zero 
            if Vr(end,j)~=0
                qreversal=qv(j);
                qfound=1;
            end
        end
    end
end

if ~isnan(qreversal)
    Ureversal=sqrt(2*qreversal/rho);
    fprintf('The first control reversal speed is %f \n',Ureversal)
else
    fprintf('No control reversal speed found\n')
end

