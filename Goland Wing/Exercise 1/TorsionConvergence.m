%This script performs a convergence analysis of the Nd first uncoupled 
% torsion modes of the Goland wing. It compares two things: 
% 1) Compares Natural Frequencies with Exact Results: uses relative error. 
% 2) Compares mode shapes with exact mode shapes, quantified with L2 norm. 

addpath("FEM\")

%% Define number of elements to sweep and desired number of modes: 
Nd=6;
Nev=5:20;

%% Prepare exact values for natural frequencies: 
L=6.096;        %[m]
GJ=987600;      %[N m^2]
Itheta=7.452;   %[kg m]  Take it just like this for this problem, for the coupled one this is Io and needs to be augmented with ms^2.

%Frequencies taken from professor's material
Betas=(1:2:11)'*pi/2/L;

Exact_Frequencies=Betas*sqrt(GJ/Itheta);

RelFreqError=nan(Nd,length(Nev));
L2error=nan(Nd,length(Nev));
for i=1:length(Nev)
    [fr,qr,FEMdata]=FEM_UncoupledTorsion(Nev(i),Nd);
    RelFreqError(:,i)=abs(fr-Exact_Frequencies)./Exact_Frequencies;
    for j=1:Nd
        SquareError=@(y) (EvaluateTorsionModalDisplacement(y,qr(:,j),FEMdata.deltay,FEMdata.Connect)/EvaluateTorsionModalDisplacement(L,qr(:,j),FEMdata.deltay,FEMdata.Connect)-EvaluateExactTorsionMode(y,Betas(j),L)).^2;
        L2error(j,i)=sqrt(1/L*integral(SquareError,0,L)); %Introduce 1/L to make this non-dimensional
    end
end


%Plot evolution of natural frequency error:
figure(1)
for j=1:Nd
    curvelabel="Mode" + j;
    plot(Nev,100*RelFreqError(j,:),'DisplayName',curvelabel)
    hold on
end
plot(Nev,ones(size(Nev)),'k--','HandleVisibility','off')
xlabel("$N_e$",'Interpreter','latex')
ylabel("Frequency Relative Error $[\%]$",'Interpreter','latex')
grid minor
legend()

%Plot evolution of modal shape error: 
figure(2)
for j=1:Nd
    curvelabel="Mode" + j;
    plot(Nev,100*L2error(j,:),'DisplayName',curvelabel)
    hold on
end
plot(Nev,ones(size(Nev)),'k--','HandleVisibility','off')
xlabel("$N_e$",'Interpreter','latex')
ylabel("Modal Shape Error $[\%]$",'Interpreter','latex')
grid minor
legend()


