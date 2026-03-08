function [fr,qr,FEMdata]=FEM_CoupledBendingTorsion(Ne,Nd)
%This function wrapper uses the procedure defined in CoupledBendingTorsion.m 
%to compute coupled bending-torsion natural frequencies and modes. It receives
%a prescribed number of elements (which are taken as joint bending-torsion
%elements) and the number of modes of interest. 
L=6.096;        %[m]
c=1.829;        %[m]
shat=0.1;          %[-] Distance between elastic and mass axes, referred to the chord. 
EI=9.77*1e6;    %[N*m2]
m=35.72;        %[kg/m]
GJ=987600;      %[N m^2]
Io=7.452;       %[kg m] 
s=shat*c;
Itheta=Io+m*s^2;

%Number of bending nodes and dofs:
NnodesB=Ne+1;
NdofsB=2*NnodesB;
%Number of torsion nodes and dofs:
NnodesT=2*Ne+1;
NdofsT=NnodesT;
%Overall dofs: 
Ndofs=NdofsB+NdofsT;

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

%Append a null rows to account for null degrees of freedom

%Reintroduce the eliminated bending dofs first:
V=[zeros(2,Ndofs-3);Vr];

%Now reintroduce the eliminated torsion displacement: 
V=[V(1:NdofsB,:);zeros(1,Ndofs-3);V(NdofsB+1:end,:)];

%Prepare returns: 

fr=Frequencies(1:Nd);
%Return eigenvectors normalized by mass matrix: 
qr=nan(Ndofs,Nd);
for i=1:Nd
    qr(:,i)=V(:,i)/sqrt(V(:,i)'*MM*V(:,i));
end

FEMdata.deltay=deltay;
FEMdata.Connect=Connect; 

end