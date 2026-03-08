function w=EvaluateTorsionModalDisplacement(y,q,deltay,Connect)
%This function evaluates the displacement at point y associated to a mode
%given by the coefficients in vector q. It needs to know deltay and the
%connectivity matrix to know the element the point is in. It uses shape
%functions for interpolation, as is conceptually consistent with the finite
%element method. Note at nodes it doesn't matter which element is
%considered, since they will be directly related to one element of q. I've
%prepared it to be able to handle a vector input x for latter applications.

%Note: q must be a column vector!



[Ne,~]=size(Connect);
w=nan(size(y));
for k=1:length(y)
    %Identify element
    EID=floor(y(k)/deltay)+1; %EID is the ID of the element where x is located. 
    
    if EID==Ne+1 %to catch fringe case where y=L
        EID=EID-1;
    end

    %Extract relevant degrees of freedom: 
    NIDs=Connect(EID,:); %ID's of the nodes relevant to this element, in the ordering I've followed so far. 

    %Evaluate shape functions:
    xi=(y(k)-(EID-1)*deltay)/deltay; 

    psi1=1-3*xi+2*xi^2;
    psi2=4*xi-4*xi^2;
    psi3=2*xi^2-xi;


    w(k)=[psi1,psi2,psi3]*q(NIDs); %to account for the change of variables we did in the rotational dofs.
end
end