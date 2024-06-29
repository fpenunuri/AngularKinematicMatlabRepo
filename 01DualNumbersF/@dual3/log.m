%LOG    Natural logarithm.
% The elements of X can be dual3 numbers
function fr = log(g)
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  
  fr0 = log(g0);
  fr1 = g1./g0;
  fr2 = g2./g0 - (g1./g0).^2;
  fr3 = (2.*g1.^3)./g0.^3 - (3.*g1.*g2)./g0.^2 + g3./g0;
  
  fr = dual3(fr0,fr1,fr2,fr3);
end
