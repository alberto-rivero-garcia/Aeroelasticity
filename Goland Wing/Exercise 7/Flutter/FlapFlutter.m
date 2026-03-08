%This script performs the flutter analysis for the Goland wing with flap.
%It returns the value of UF and the frequency versus damping plots for a
%given number of elements and size of the modal basis. 
addpath("Auxiliary Functions\")
Ne=15; 
Nd=6; 
Umax=200; 
Nu=2000;

UF=FlapPKGoland(Ne,Nd,Umax,Nu,1);

display(UF)