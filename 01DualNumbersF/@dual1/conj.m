%CONJ   Complex conjugate.
%   CONJ(X) is the complex conjugate of X.
%   For a complex X, CONJ(X) = REAL(X) - i*IMAG(X).
% The elements of X can be dual1 numbers
% X = f0 + f1 e1, 
% conj(X) = conj(f0) + conj(f1) e1 
function fr = conj(g) 
  g0 = g.f0;
  g1 = g.f1;
  
  fr0 = conj(g0);
  fr1 = conj(g1);
  
  fr = dual1(fr0,fr1);
end

