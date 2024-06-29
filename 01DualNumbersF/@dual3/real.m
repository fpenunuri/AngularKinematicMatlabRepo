%   REAL(X) is the real part of X.
% dual3 case
function fr = real(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  
  fr0 = real(g0);
  fr1 = real(g1);
  fr2 = real(g2);
  fr3 = real(g3);
  
  fr = dual3(fr0,fr1,fr2,fr3);
end