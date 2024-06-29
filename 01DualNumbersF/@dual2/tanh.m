%TANH   Hyperbolic tangent.
% The elements of X can be dual2 numbers
function fr = tanh(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  
  fr0 = tanh(g0);
  fr1 = g1./cosh(g0).^2;
  fr2 = g2./cosh(g0).^2 - (2.*g1.^2).*tanh(g0)./cosh(g0).^2;
  
  fr = dual2(fr0,fr1,fr2);
end
 