%CONJ   Complex conjugate.
%   CONJ(X) is the complex conjugate of X.
%   For a complex X, CONJ(X) = REAL(X) - i*IMAG(X).
% The elements of X can be dual3 numbers
% X = f0 + f1 e1 + f2 e2 + f3 e3, 
% conj(X) = conj(f0) + conj(f1) e1 + conj(f2) e2 + conj(f3) e3
function fr = conj(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  
  fr0 = conj(g0);
  fr1 = conj(g1);
  fr2 = conj(g2);
  fr3 = conj(g3);
  
  fr = dual3(fr0,fr1,fr2,fr3);
end

