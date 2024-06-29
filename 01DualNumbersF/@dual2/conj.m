%CONJ   Complex conjugate.
%   CONJ(X) is the complex conjugate of X.
%   For a complex X, CONJ(X) = REAL(X) - i*IMAG(X).
% The elements of X can be dual2 numbers
% X = f0 + f1 e1 + f2 e2, 
% conj(X) = conj(f0) + conj(f1) e1 + conj(f2) e2
function fr = conj(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  
  fr0 = conj(g0);
  fr1 = conj(g1);
  fr2 = conj(g2);
  
  fr = dual2(fr0,fr1,fr2);
end

