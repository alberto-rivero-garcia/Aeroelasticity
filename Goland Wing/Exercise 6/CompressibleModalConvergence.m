%This script studies the evolution of the flutter velocity prediction
%obtained with the PK method for the Goland wing depending on the size of
%the modal basis. For plotting, see the PlotCompressiblePKGoland.m
%auxiliary function. 
addpath("Auxiliary Functions\")
%Set up basic parameters: 
Ne=10;
Ndv=2:6;
Umax=200;
Nu=1000;

%Initialize UF vector and perform flutter iteration. 
UFv=nan(size(Ndv));
UFv2=nan(size(Ndv));

for i=1:length(Ndv)
    UFv(i)=CompressiblePKGoland(Ne,Ndv(i),Umax,Nu);
    UFv2(i)=CompressiblePKGoland(Ne,2*Ndv(i),Umax,Nu);
end

%To study modal convergence, I will actually employ Richardson's
%extrapolation technique. (assuming first order convergence). 

UFRichardson=2*UFv2-UFv;

Epsilon=abs(UFv-UFRichardson)./UFRichardson; %we can consider it converged even with 2 modes, it gets better afterwards though. 

