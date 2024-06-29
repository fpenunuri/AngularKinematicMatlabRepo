%LOG    Natural logarithm.
% The elements of X can be dual2 numbers
function fr = log(g)
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  
  fr0 = log(g0);
  fr1 = g1./g0;
  fr2 = g2./g0 - (g1./g0).^2;
  
  fr = dual2(fr0,fr1,fr2);
end