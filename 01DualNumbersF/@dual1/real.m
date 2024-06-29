%   REAL(X) is the real part of X.
%
function fr = real(g) 
  g0 = g.f0;
  g1 = g.f1;
  
  fr0 = real(g0);
  fr1 = real(g1);
  
  fr = dual1(fr0,fr1);
end