%Check dCkdk: based on fundamental theorem of calculus: 
%All checks passed!
k_v=linspace(0.1,1,100);
Ck_v=Ck(k_v);

%Check 1: real part: 
realCk_v_int=nan(size(k_v));
for i=1:length(k_v)
    realCk_v_int(i)=real(Ck_v(1))+integral(@(k)real(dCkdk(k)),k_v(1),k_v(i));
end

figure(1)
plot(k_v,real(Ck_v),'DisplayName','Exact')
hold on
plot(k_v,realCk_v_int,'DisplayName','Integrated')
legend()

%Check 2: imaginary part: 
imagCk_v_int=nan(size(k_v));
for i=1:length(k_v)
    imagCk_v_int(i)=imag(Ck_v(1))+integral(@(k)imag(dCkdk(k)),k_v(1),k_v(i));
end

figure(2)
plot(k_v,imag(Ck_v),'DisplayName','Exact')
hold on
plot(k_v,imagCk_v_int,'DisplayName','Integrated')
legend()


%Prepare auxilary functions: 
function C=Ck(k)
%This function computes Theodorsen's function from a general vector input
%k. Returns a complex vector. 
C=nan(size(k));

for i=1:length(k)
    C(i)=besselh(1,2,k(i))/(besselh(1,2,k(i))+1i*besselh(0,2,k(i)));
end
end

function dCk=dCkdk(k)
%This function computes the derivative of Theodorsen's function wrt k,
%returning it as a complex number. Takes function input to admit
%integration with matlab's integral command. 
dCk=nan(size(k));

for i=1:length(k)
    dH1dk=1/2*(besselh(0,2,k(i))-besselh(2,2,k(i)));
    dH0dk=1/2*(besselh(-1,2,k(i))-besselh(1,2,k(i)));
    dCk(i)=1i*(dH1dk*besselh(0,2,k(i))-besselh(1,2,k(i))*dH0dk)/((besselh(1,2,k(i))+1i*besselh(0,2,k(i)))^2);
end
end
