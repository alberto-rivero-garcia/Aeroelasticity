function w=EvaluateExactBendingMode(y,beta,L)
%This function evaluates the vertical displacement of the exact bending
%mode associated with the input beta, for a beam of length L. It is
%prepared to provide a unit displacement on the free edge of the beam. This
%is not a useful choice in terms of further implementation, but it is very 
%useful to contrast it with the modes. Capable of handling vector inputs
%for implementation into other workflows. Will return a vector that matches
%the dimensions of the input. (very relevant for that last point). 
A1=(cosh(beta*L)-cos(beta*L)+(cos(beta*L)+cosh(beta*L))/(sin(beta*L)+sinh(beta*L))*(sin(beta*L)-sinh(beta*L)))^(-1);

w=A1*(cosh(beta*y)-cos(beta*y)+(cos(beta*L)+cosh(beta*L))/(sin(beta*L)+sinh(beta*L))*(sin(beta*y)-sinh(beta*y)));
end