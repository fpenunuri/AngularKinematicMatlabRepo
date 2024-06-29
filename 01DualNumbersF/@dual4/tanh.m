%TANH   Hyperbolic tangent.
% The elements of X can be dual4 numbers
function fr = tanh(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  g4 = g.f4;
  
  vd3 = tanh(dual3(g0,g1,g2,g3));  
  vd4 = (g4 + 12.*g1.^2 .* g2.*(2 - 3./cosh(g0).^2) + 2.*(-3.*g2.^2 - ...
         4.*g1.*(g1.^3 + g3) + (12.*g1.^4)./cosh(g0).^2).*tanh(g0))./ ...
         cosh(g0).^2;
  
  fr = dual4(vd3.f0,vd3.f1,vd3.f2,vd3.f3,vd4);
end
