function Ka_local=TorsionLocalAeroMatrix(ymin,ymax,Kast)
%This function performs the integration necessary to construct the local
%aerodynamic stiffness matrix for control reversal analysis:
%Define xi:
deltay=ymax-ymin;
xi=@(y)(y-ymin)/deltay;

%Define all shape functions: 
%Torsion shape functions: 
psit1=@(y)1-3*xi(y)+2*xi(y).^2;
psit2=@(y)4*xi(y)-4*xi(y).^2;
psit3=@(y)2*xi(y).^2-xi(y);


%Start integrating: 
Ka_tt=integral(@(y)[psit1(y);psit2(y);psit3(y)]*Kast(1,1)*[psit1(y),psit2(y),psit3(y)],ymin,ymax,'ArrayValued',true);
Ka_tb=integral(@(y)[psit1(y);psit2(y);psit3(y)]*Kast(1,2)*Betafun(y,1),ymin,ymax,'ArrayValued',true);
Ka_bt=integral(@(y)Betafun(y,1)*Kast(2,1)*[psit1(y),psit2(y),psit3(y)],ymin,ymax,'ArrayValued',true);
Ka_bb=integral(@(y)Betafun(y,1)*Kast(2,2)*Betafun(y,1),ymin,ymax,'ArrayValued',true);

%Construct local aero matrix: 
Ka_local=[Ka_tt,Ka_tb;Ka_bt,Ka_bb];

end