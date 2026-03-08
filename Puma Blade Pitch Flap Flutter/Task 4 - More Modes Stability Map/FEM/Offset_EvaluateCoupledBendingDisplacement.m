function w=Offset_EvaluateCoupledBendingDisplacement(x,q,deltax,Connect,xmin)
%This function evaluates the displacement at point y associated to a mode
%given by the coefficients in vector q. It needs to know deltay and the
%connectivity matrix to know the element the point is in. It uses shape
%functions for interpolation, as is conceptually consistent with the finite
%element method. Note at nodes it doesn't matter which element is
%considered, since they will be directly related to one element of q. I've
%prepared it to be able to handle a vector input x for latter applications.

%Note: q must be a column vector!

%Needs to be corrected to properly account for the fact that elements start
%in x=xe.

[Ne,~]=size(Connect);
w=nan(size(x));
for k=1:length(x)
    %Identify element
    EID=floor((x(k)-xmin)/deltax)+1; %EID is the ID of the element where x is located. 
    
    if EID==Ne+1 %to catch fringe case where y=L
        EID=EID-1;
    end

    %Extract relevant degrees of freedom: 
    NIDs=Connect(EID,1:4); %ID's of the nodes relevant to this element, in the ordering I've followed so far. 

    %Evaluate shape functions:
    xi=(x(k)-xmin-(EID-1)*deltax)/deltax; 

    psi1=1-3*xi^2+2*xi^3;
    psi2=xi-2*xi^2+xi^3;
    psi3=3*xi^2-2*xi^3;
    psi4=-xi^2+xi^3;

    w(k)=[psi1,psi2*deltax,psi3,psi4*deltax]*q(NIDs); %to account for the change of variables we did in the rotational dofs.
end
end