%COSH   Hyperbolic cosine.
% The elements of X can be dual2 numbers
function fr = cosh(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  
  fr0 = cosh(g0);
  fr1 = sinh(g0).*g1;
  fr2 = g2.*sinh(g0) + cosh(g0).*(g1.^2);
  
  fr = dual2(fr0,fr1,fr2);
end
