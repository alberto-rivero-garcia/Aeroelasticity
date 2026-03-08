function NA=Free_Play_SIDF(A,delta,k)
%This function implements the sinusoidal input describing function for a
%freeplay of range 2delta and stiffness k. 
NA=nan(size(A));

for i=1:length(A)
    if A(i)<=delta
        NA(i)=0;
    else
        NA(i)=k*(1-2/pi*(asin(delta/A(i))+delta/A(i)*sqrt(1-(delta/A(i))^2)));
    end
end
end