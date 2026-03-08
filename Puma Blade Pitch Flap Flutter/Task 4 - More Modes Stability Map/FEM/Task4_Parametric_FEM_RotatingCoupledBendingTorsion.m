function [fr,qr,FEMdata,FEMCheck]=Task4_Parametric_FEM_RotatingCoupledBendingTorsion(Ne,Nd,Ktheta,f_eps,eps_waypoints)
%This function returns the Nd first modal frequencies and mode shapes of 
% the rotating Puma blade, using Ne elements in the discretization. It also
% returns relevant FEMdata for results plotting. 

%It has been developed from the v2 version of the original FEM, and is 
%prepared to take control chain stiffness and desired position of the blade
%cog (in the chord frame) as inputs.  

%v3 of the parametric implementation also adjusts the torsional inertia of
%the section, which leads to a PD mass matrix. For that, it leverages v2 of
%the RotatingElement function. Furthermore, it includes a sanity check for
%eigenvalues that is used to discard cases where the system is unstable in
%vacuo. 

xe= 0.289;      %[m]
R=7.490;        %[m]

Omega=4.5*2*pi; %[rad/s] From Quaranta's code, not specified in the data documentation




%Number of bending nodes and dofs:
NnodesB=Ne+1;
NdofsB=2*NnodesB;
%Number of torsion nodes and dofs:
NnodesT=2*Ne+1;
NdofsT=NnodesT;
%Overall dofs: 
Ndofs=NdofsB+NdofsT;

%Element size: 
deltax=(R-xe)/Ne;

%Define position of element boundaries:  
b_pos=xe:deltax:R;

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

%Assemble global mass and stiffness matrices. 
KK=zeros(Ndofs,Ndofs);
MM=zeros(Ndofs,Ndofs);
for i=1:Ne
    %Construct local element matrices:  
    [KKlocal,MMlocal]=Task4_Parametric_RotatingLocalElementMatrices(b_pos(i),b_pos(i+1),Omega,R,f_eps,eps_waypoints);
    %Add contributions to the global matrices: 
    KK(Connect(i,:),Connect(i,:))=KK(Connect(i,:),Connect(i,:))+KKlocal;
    MM(Connect(i,:),Connect(i,:))=MM(Connect(i,:),Connect(i,:))+MMlocal;    
end

%Add term related to the Robin BC at the root. Applied at the first torsion
%dof of the first element, i.e, torsion at x=e. 
KK(Connect(1,5),Connect(1,5))=KK(Connect(1,5),Connect(1,5))+Ktheta;

%Eliminate vertical displacement of the root (easy to see
%because it is the first dof!)
MMr=MM(2:end,2:end);
KKr=KK(2:end,2:end);
%I've checked both of this are full-rank -> so there's no singularity
%issues here. 

%Solve eigenvalue-eigenvector problem: this requires some reframing, since
%it is returning lambda^2 and it's much more difficult to choose the
%appropriate eigenvalues here. It seems that a different approach might be
%worth it. 
[Vr,D]=eigs(-KKr,MMr,Nd,'smallestabs'); 


[Frequencies,sorting]=sort(imag(sqrt(diag(D))));

Vr=Vr(:,sorting); %some bona fide matlab magic to order the matrix by columns according to the sorting that came out of the previous step

%Append null row to account for null degree of freedom
V=[zeros(1,Nd);Vr];

%Prepare returns: 
fr=Frequencies(1:Nd);
%Return eigenvectors normalized by mass matrix: 
qr=nan(Ndofs,Nd);
for i=1:Nd
    qr(:,i)=V(:,i)/sqrt(V(:,i)'*MM*V(:,i)); 
end

FEMdata.deltax=deltax;
FEMdata.Connect=Connect; 
FEMdata.xmin=b_pos(1);
FEMdata.xmax=b_pos(end);

%Sanity check on eigenvalues: 
raw_eigs=sqrt(diag(D));
sorted_eigs=raw_eigs(sorting);
if all(real(sorted_eigs(1:Nd))==0) && all(imag(sorted_eigs(1:Nd))>0)
    FEMCheck=1;
else
    FEMCheck=0;
end

end