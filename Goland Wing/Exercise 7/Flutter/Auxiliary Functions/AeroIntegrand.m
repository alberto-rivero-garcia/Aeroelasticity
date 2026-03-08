function val=AeroIntegrand(y,MyModei,MyModej,matrix)
%This function is prepared to act as an integrand for the computation of
%the aerodynamic modal projection integrals. This wrapper is needed to make it able to handle vector
%inputs and thus, use the integral command. 

val=nan(size(y));
for i=1:length(y)
    val(i)=MyModei(y(i))'*matrix*[-1,0,0;0,1,0;0,0,1]*MyModej(y(i)); %change here to account for the negative sign in psiwj
end
end