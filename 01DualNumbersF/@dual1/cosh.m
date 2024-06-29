%COSH   Hyperbolic cosine.
% The elements of X can be dual1 numbers
function fr = cosh(g) 
  g0 = g.f0;
  g1 = g.f1;
  
  fr0 = cosh(g0);
  fr1 = sinh(g0).*g1;
  
  fr = dual1(fr0,fr1);
end
