%This script estimates the value of ko in the Richardson extrapolation
%technique (see wikipedia and paper notebook), for the first six coupled
%natural frequencies. It does so starting at a specified number of
%elements, to see how it affects the order of convergence. It will be the
%basis for determining the value of k0 to be used in Richardson
%extrapolation for exercise 2. It uses wikipedia's method for ko
%estimation, with t=2 and thus s=4. 

%Select baseline number of elements: 
Ne=10; 
Nd=6; 

%Evaluate problem at different values of h.
[Ah,~,~]=FEM_CoupledBendingTorsion(Ne,Nd);
[Ah_t,~,~]=FEM_CoupledBendingTorsion(2*Ne,Nd);
[Ah_t2,~,~]=FEM_CoupledBendingTorsion(4*Ne,Nd);

%Compute estimates for tko:
tkov=(Ah-Ah_t)./(Ah_t-Ah_t2);

%Return ko estimations: 
ko=log2(tkov)
