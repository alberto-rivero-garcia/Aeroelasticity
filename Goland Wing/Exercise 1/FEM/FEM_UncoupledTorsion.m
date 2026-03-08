function [fr,qr,FEMdata]=FEM_UncoupledTorsion(Ne,Nd)
%This function wrapper uses the procedure defined in UncoupledTorsion.m to
%compute uncoupled bending natural frequencies and modes. It then returns
%those to be compared with analytical solutions and evaluate convergence. 

L=6.096;        %[m]
GJ=987600;      %[N m^2]
Itheta=7.452;   %[kg m]  Take it just like this for this problem, for the coupled one this is Io and needs to be augmented with ms^2. Explain it in the report!!!!!

Nnodes=2*Ne+1; %First element provides 3 nodes and from then on it's 2 more per element. 
%Define element size and node locations: 
deltay=L/Ne;

%Define connectivity matrix: 
Ndofs=Nnodes; 
Connect=nan(Ne,3);

for i=1:Ne
    Connect(i,:)=[1+2*(i-1),2+2*(i-1),3+2*(i-1)];
end

K_local=GJ/deltay/3*[7,-8,1;-8,16,-8;1,-8,7];
M_local=Itheta*deltay/30*[4,2,-1;2,16,2;-1,2,4];


KK=zeros(Ndofs,Ndofs);
MM=zeros(Ndofs,Ndofs);

for i=1:Ne
    KK(Connect(i,:),Connect(i,:))=KK(Connect(i,:),Connect(i,:))+K_local;
    MM(Connect(i,:),Connect(i,:))=MM(Connect(i,:),Connect(i,:))+M_local;
end

MMr=MM(2:end,2:end);
KKr=KK(2:end,2:end);

[Vr,D]=eig(-KKr,MMr);

[Frequencies,sorting]=sort(imag(sqrt(diag(D))));

Vr=Vr(:,sorting); %some bona fide matlab magic to order the matrix by columns according to the sorting that came out of the previous step

%Append a null rows at the top of V to account for the restricted dofs
%(they are still relevant for plotting)

V=[zeros(1,Ndofs-1);Vr];

%Return desired frequencies and mode shapes, as well as a FEM data structure:
fr=Frequencies(1:Nd);
qr=V(:,1:Nd);
FEMdata.deltay=deltay;
FEMdata.Connect=Connect; 

end
