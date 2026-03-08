function val=Betafun(y,q)
%This function returns the section's beta value given the corresponding
%mode. 
yL=3.35; %[m]
yU=3.35+1.8; %[m]

val=nan(size(y));

for i=1:length(y)
    if yL<=y(i) && y(i)<=yU
        val(i)=q(end);
        else
        val(i)=0;
    end
end

end