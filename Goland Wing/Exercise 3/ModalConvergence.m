%This script studies the evolution of the flutter velocity prediction
%obtained with the PK method for the Goland wing depending on the size of
%the modal basis. For plotting of frequency and damping versus airspeed
%plots, look at the plotting auxiliary function. 
addpath("Auxiliary Functions\")

%Set up basic parameters: 
Ne=15;
Ndv=2:6;
Umax=200;
Nu=1000;

%Initialize UF vector and perform flutter iteration. 
UFv=nan(size(Ndv));

for i=1:length(Ndv)
    UFv(i)=PKGoland(Ne,Ndv(i),Umax,Nu);
end

%From Goland's paper we see that the exact flutter speed is reported as 307
%mph. This is, approximately:
UFexact=137.241; %[m/s]

%Let's do a plot for the relative error as a function of the number of
%modes: 
Epsilon=abs(UFv-UFexact)/UFexact;

figure(1)
plot(Ndv,100*Epsilon)
hold on
grid minor
%plot(Ndv,ones(size(Ndv)),'r--')
xlabel('Number of Modes $[-]$','Interpreter','latex')
ylabel('Relative Error In $U_F \; [\%]$','Interpreter','latex')

%We see that two modes are more than enough, even though the solution later
%stabilizes to something a bit higher. In any case, this is perfect. 