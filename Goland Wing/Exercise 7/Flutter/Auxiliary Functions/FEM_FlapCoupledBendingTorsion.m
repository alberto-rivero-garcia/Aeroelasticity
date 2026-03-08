function [fr,qr,FEMdata]=FEM_FlapCoupledBendingTorsion(Ne,Nd)
%This function wrapper uses a procedure derived from CoupledBendingTorsion.m 
%to compute coupled bending-torsion natural frequencies and modes. It receives
%a prescribed number of elements (which are taken as joint bending-torsion
%elements) and the number of modes of interest. It also includes the
%effects of the flap, in preparation for a more advanced flutter analysis. 

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


%Number of bending nodes and dofs:
NnodesB=Ne+1;
NdofsB=2*NnodesB;
%Number of torsion nodes and dofs:
NnodesT=2*Ne+1;
NdofsT=NnodesT;
%Overall dofs: 
Ndofs=NdofsB+NdofsT+1; %now we must also consider Beta as a dof.

%Element size: 
deltay=L/Ne;

%Define connectivity matrix: 
%Now we link each element to its degrees of freedom.
%Each element has 7 dofs. The first 4 are bending, while the next 3 are
%torsion. I will order the generalized coordinates by putting all bending
%dofs first and then putting the torsion ones. 

Connect=nan(Ne,7);

for i=1:Ne
    %Assign Bending dofs: see the uncoupled script for an explanation
    Connect(i,1:4)=[1+(i-1)*2,2+(i-1)*2,3+(i-1)*2,4+(i-1)*2];
    %Assign Torsion dofs: the key is in the shift by NdofsB
    Connect(i,5:7)=NdofsB+[1+2*(i-1),2+2*(i-1),3+2*(i-1)];
end

%Bending Problem:
Kww_local=EI/deltay^3*[12,6*deltay,-12,6*deltay;6*deltay,4*deltay^2,-6*deltay,2*deltay^2;...
    -12,-6*deltay,12,-6*deltay;6*deltay,2*deltay^2,-6*deltay,4*deltay^2];


Mww_local=m*deltay/420*[156,22*deltay,54,-13*deltay;22*deltay,4*deltay^2,13*deltay,-3*deltay^2;...
    54,13*deltay,156,-22*deltay;-13*deltay,-3*deltay^2,-22*deltay,4*deltay^2];

%Torsion problem:
Ktt_local=GJ/deltay/3*[7,-8,1;-8,16,-8;1,-8,7];
Mtt_local=Itheta*deltay/30*[4,2,-1;2,16,2;-1,2,4];

%Mass coupling: 
Mwt_local=-m*deltay*s/60*[11,20,-1;deltay,4*deltay,0;-1,20,11;0,-4*deltay,-deltay];

%Coupled Bending-Torsion Problem local element: 
K_local=[Kww_local,zeros(4,3);zeros(3,4),Ktt_local];
M_local=[Mww_local,Mwt_local;Mwt_local',Mtt_local]; %note the matrix will still be symmetric -see the transpose on the 2-1 position

KK=zeros(Ndofs,Ndofs);
MM=zeros(Ndofs,Ndofs);

for i=1:Ne
    KK(Connect(i,:),Connect(i,:))=KK(Connect(i,:),Connect(i,:))+K_local;
    MM(Connect(i,:),Connect(i,:))=MM(Connect(i,:),Connect(i,:))+M_local;
end

%Modify mass matrix to account for flap: (see notebook for model
%explanation and development).
%Identify yf's element: 
EID=floor(yf/deltay)+1; %EID is the ID of the element where x is located. 
NIDs=Connect(EID,:);    %Extract nodal IDs

%Evaluate shape functions: 
%Bending shape functions: 
xi=(yf-(EID-1)*deltay)/deltay; 
psiw1=1-3*xi^2+2*xi^3;
psiw2=(xi-2*xi^2+xi^3)*deltay;
psiw3=3*xi^2-2*xi^3;
psiw4=(-xi^2+xi^3)*deltay;
psiw=[psiw1,psiw2,psiw3,psiw4]';
%Torsion shape functions: 
psit1=1-3*xi+2*xi^2;
psit2=4*xi-4*xi^2;
psit3=2*xi^2-xi;
psit=[psit1,psit2,psit3]';

%Add contributions: 
%Relating bending and torsion
MM(NIDs(1:4),NIDs(1:4))=MM(NIDs(1:4),NIDs(1:4))+psiw*Mf*psiw';
MM(NIDs(5:7),NIDs(5:7))=MM(NIDs(5:7),NIDs(5:7))+psit*(If+(xf-xe)^2*Mf)*psit';
MM(NIDs(1:4),NIDs(5:7))=MM(NIDs(1:4),NIDs(5:7))-psiw*Mf*(xf-xe)*psit';
MM(NIDs(5:7),NIDs(1:4))=MM(NIDs(5:7),NIDs(1:4))-psit*Mf*(xf-xe)*psiw';
%Coupling bending and torsion with flap deflection: 
MM(NIDs(1:4),end)=MM(NIDs(1:4),end)-psiw*Mf*cf/2;
MM(NIDs(5:7),end)=MM(NIDs(5:7),end)+psit*(If+Mf*cf/2*(xf-xe));

%Introduce the beta equation: 
%Direct mass and stiffness
MM(end,end)=Ih;
KK(end,end)=KBeta;
%Coupling terms: 
MM(end,NIDs(1:4))=MM(end,NIDs(1:4))-Mf*cf/2*psiw';
MM(end,NIDs(5:7))=MM(end,NIDs(5:7))+(If+Mf*cf/2*(xf-xe))*psit';

%Apply BCs

%Eliminate torsion of the root: (I do it first because it's easier to find
%like this). Eliminate the first torsional dof. 
MMr=MM([1:NdofsB, NdofsB+2:end],[1:NdofsB, NdofsB+2:end]);
KKr=KK([1:NdofsB, NdofsB+2:end],[1:NdofsB, NdofsB+2:end]);

%Eliminate vertical displacement and bending angle of the root (easy to see
%because they are the two first dofs)
MMr=MMr(3:end,3:end);
KKr=KKr(3:end,3:end);

%Solve eigenvalue-eigenvector problem: 
[Vr,D]=eig(-KKr,MMr);

[Frequencies,sorting]=sort(imag(sqrt(diag(D))));

Vr=Vr(:,sorting); %some bona fide matlab magic to order the matrix by columns according to the sorting that came out of the previous step

%Append null rows to account for null degrees of freedom

%Reintroduce the eliminated bending dofs first:
V=[zeros(2,Ndofs-3);Vr];

%Now reintroduce the eliminated torsion displacement: 
V=[V(1:NdofsB,:);zeros(1,Ndofs-3);V(NdofsB+1:end,:)];

%Prepare returns: 

fr=Frequencies(1:Nd);
%Return eigenvectors normalized by mass matrix: 
qr=nan(Ndofs,Nd);
for i=1:Nd
    qr(:,i)=V(:,i)/sqrt(V(:,i)'*MM*V(:,i)); %There's an issue here!!!!!
end

FEMdata.deltay=deltay;
FEMdata.Connect=Connect; 

end