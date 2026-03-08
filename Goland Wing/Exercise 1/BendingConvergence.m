%This script performs a convergence analysis of the Nd first uncoupled 
% bending modes of the Goland wing. It compares two things: 
% 1) Compares Natural Frequencies with Exact Results: uses relative error. 
% 2) Compares mode shapes with exact mode shapes, quantified with L2 norm. 
addpath("FEM\")
%% Define number of elements to sweep and desired number of modes: 
Nd=6;
Nev=5:20;

%% Prepare exact values for natural frequencies: 
L=6.096; %[m]
EI=9.77*1e6; %[N*m2]
m=35.72; %[kg/m]
Betas=[0.596864,1.49418,2.50025,3.49999,4.5,5.5]'*pi/L; %taken from professor's material 
Exact_Frequencies=Betas.^2*sqrt(EI/m);

RelFreqError=nan(Nd,length(Nev));
L2error=nan(Nd,length(Nev));
for i=1:length(Nev)
    [fr,qr,FEMdata]=FEM_UncoupledBending(Nev(i),Nd);
    RelFreqError(:,i)=abs(fr-Exact_Frequencies)./Exact_Frequencies;
    for j=1:Nd
        SquareError=@(y) (EvaluateBendingModalDisplacement(y,qr(:,j),FEMdata.deltay,FEMdata.Connect)/EvaluateBendingModalDisplacement(L,qr(:,j),FEMdata.deltay,FEMdata.Connect)-EvaluateExactBendingMode(y,Betas(j),L)).^2;
        L2error(j,i)=sqrt(1/L*integral(SquareError,0,L)); %Introduce 1/L to make this non-dimensional
    end
end


%Plot evolution of natural frequency error: note natural frequencies are
%returned in rad/s!
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


 


