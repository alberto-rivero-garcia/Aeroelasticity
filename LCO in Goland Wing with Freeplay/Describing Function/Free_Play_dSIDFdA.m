function dNA=Free_Play_dSIDFdA(A,delta,k)
%This function implements the derivative wrt A of the sinusoidal input describing
%function for a free play of range 2delta and stiffness k. 

dNA=nan(size(A));

for i=1:length(A)
    if A(i)<=delta
        dNA(i)=0;
    else
        dNA(i)=4*k*delta/pi/(A(i)^2)*sqrt(1-(delta/A(i))^2);
    end
end
end