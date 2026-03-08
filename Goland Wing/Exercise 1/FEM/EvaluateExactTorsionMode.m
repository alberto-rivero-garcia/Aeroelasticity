function w=EvaluateExactTorsionMode(y,beta,L)
%This function evaluates the vertical displacement of the exact bending
%mode associated with the input beta, for a beam of length L. It is
%prepared to provide a unit displacement on the free edge of the beam. This
%is not a useful choice in terms of further implementation, but it is very 
%useful to contrast it with the modes. Capable of handling vector inputs
%for implementation into other workflows. Will return a vector that matches
%the dimensions of the input. (very relevant for that last point). 
A1=1/sin(beta*L);

w=A1*sin(beta*y);
end