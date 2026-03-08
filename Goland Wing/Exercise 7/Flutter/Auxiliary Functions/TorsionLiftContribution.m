function L_local=TorsionLiftContribution(ymin,ymax,Kast)
%This function computes the contribution of a given element to the overall
%lift on the Goland wing. It assumes that the lift is represented on the
%first row of the static aerodynamic matrix. 
%Define xi:
deltay=ymax-ymin;
xi=@(y)(y-ymin)/deltay;

%Define all shape functions: 
%Torsion shape functions: 
psit1=@(y)1-3*xi(y)+2*xi(y).^2;
psit2=@(y)4*xi(y)-4*xi(y).^2;
psit3=@(y)2*xi(y).^2-xi(y);

%Start integrating: 
L_t=integral(@(y)Kast(1,2)*[psit1(y),psit2(y),psit3(y)],ymin,ymax,'ArrayValued',true);
L_b=integral(@(y)Kast(1,3)*Betafun(y,1),ymin,ymax,'ArrayValued',true);

L_local=[L_t,L_b];
end