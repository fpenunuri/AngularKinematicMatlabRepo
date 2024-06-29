%SINH   Hyperbolic sine.
% The elements of X can be dual4 numbers
function fr = sinh(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  g4 = g.f4;
  
  vd3 = sinh(dual3(g0,g1,g2,g3));  
  vd4 = (6.*g1.^2 .* g2 + g4).*cosh(g0) + (g1.^4 + 3 *g2.^2 + ...
         4.*g1.*g3).*sinh(g0);
  
  fr = dual4(vd3.f0,vd3.f1,vd3.f2,vd3.f3,vd4);
end