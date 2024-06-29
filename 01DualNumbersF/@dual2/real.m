%   REAL(X) is the real part of X.
%dual2 case
function fr = real(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  
  fr0 = real(g0);
  fr1 = real(g1);
  fr2 = real(g2);
  
  fr = dual2(fr0,fr1,fr2);
end