function [fr,qr,FEMdata]=FEM_UncoupledBending(Ne,Nd)
%This function wrapper uses the procedure defined in UncoupledBending.m to
%compute uncoupled bending natural frequencies and modes. It then returns
%those to be compared with analytical solutions and evaluate convergence. 

L=6.096; %[m]
EI=9.77*1e6; %[N*m2]
m=35.72; %[kg/m]

Nnodes=Ne+1; %First element provides 2 nodes and from then on it's one more per element. 
%Define element size and node locations: 
deltay=L/Ne;

%Define connectivity matrix: but relate it to dofs (as is consistent with the method) 
%My numbering scheme is 1 -> w1 at el. 1, 2-> phi1 at el. 1, 3-> w2 at el1,
%4-> phi2 at el. 1 and so on. 
Ndofs=2*Nnodes; %two per node
Connect=nan(Ne,4);

for i=1:Ne
    Connect(i,:)=[1+(i-1)*2,2+(i-1)*2,3+(i-1)*2,4+(i-1)*2]; %element 1 has dofs from 1 to 4. element two has dofs from 3 to 6 and so on.
end
 
K_local=EI/deltay^3*[12,6*deltay,-12,6*deltay;6*deltay,4*deltay^2,-6*deltay,2*deltay^2;...
    -12,-6*deltay,12,-6*deltay;6*deltay,2*deltay^2,-6*deltay,4*deltay^2];


M_local=m*deltay/420*[156,22*deltay,54,-13*deltay;22*deltay,4*deltay^2,13*deltay,-3*deltay^2;...
    54,13*deltay,156,-22*deltay;-13*deltay,-3*deltay^2,-22*deltay,4*deltay^2];

KK=zeros(Ndofs,Ndofs);
MM=zeros(Ndofs,Ndofs);

for i=1:Ne
    KK(Connect(i,:),Connect(i,:))=KK(Connect(i,:),Connect(i,:))+K_local;
    MM(Connect(i,:),Connect(i,:))=MM(Connect(i,:),Connect(i,:))+M_local;
end

MMr=MM(3:end,3:end);
KKr=KK(3:end,3:end);

[Vr,D]=eig(-KKr,MMr);

[Frequencies,sorting]=sort(imag(sqrt(diag(D))));

Vr=Vr(:,sorting); %some bona fide matlab magic to order the matrix by columns according to the sorting that came out of the previous step

%Append two null rows at the top of V to account for the restricted dofs
%(they are still relevant for plotting)

V=[zeros(2,Ndofs-2);Vr];

%Return desired frequencies and mode shapes, as well as a FEM data structure: 
fr=Frequencies(1:Nd);
qr=V(:,1:Nd);
FEMdata.deltay=deltay;
FEMdata.Connect=Connect; 

end