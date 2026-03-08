%This script does the basic plotting I want for the solution of exercise 2.
%These are 1) an evolution of modal natural frequencies as a function of
%Ne and 2) evolution of relative error wrt Richardson Prediction as a
%function of Ne. 

Nev=5:15;
Nd=6;
ko=3; 

fnum=nan(Nd,length(Nev));
error=nan(Nd,length(Nev));
for i=1:length(Nev)
    [fr,~,~]=FEM_CoupledBendingTorsion(Nev(i),Nd);
    [fr2,~,~]=FEM_CoupledBendingTorsion(2*Nev(i),Nd);
    fRichardson=(2^ko*fr2-fr)/(2^ko-1);
    fnum(:,i)=fr;
    error(:,i)=abs(fr-fRichardson)./fRichardson;
end


%Plot evolution of modal frequencies: 
figure(1)
for j=1:Nd
    curvelabel="Mode" + j;
    plot(Nev,fnum(j,:),'DisplayName',curvelabel)
    hold on
end
xlabel("Number of elements")
ylabel("Natural frequency $[Hz]$",'Interpreter','latex')
grid minor
legend()

%Plot Evolution of Relative error
figure(2)
for j=1:Nd
    curvelabel="Mode " + j;
    plot(Nev,100*error(j,:),'DisplayName',curvelabel)
    hold on
end
plot(Nev,ones(size(Nev)),'k--','HandleVisibility','off')
xlabel("$N_e \;[-]$",'Interpreter','latex')
ylabel("Coupled Frequency Relative Error $[\%]$",'Interpreter','latex')
grid minor
legend()


