%This script tests case c) of exercise 5 with the number of elements, modes
%and so on derived from previous studies. 
addpath("Auxiliary Functions\")
%Define parameters: 
xmhat=0.5; 
Ne=10;
Nd=6;
Umax=200;
Nu=10000;
Display=1; %set to 1 to display mode shapes and frequencies. 
%Perform PK: 
UF=PK_TipMassGoland(Ne,Nd,Umax,Nu,xmhat,Display);

display(UF)
