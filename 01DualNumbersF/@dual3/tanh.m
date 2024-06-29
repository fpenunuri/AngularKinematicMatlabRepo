%TANH   Hyperbolic tangent.
% The elements of X can be dual3 numbers
function fr = tanh(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;  
  
  fr0 = tanh(g0);
  fr1 = g1./cosh(g0).^2;
  fr2 = g2./cosh(g0).^2 - (2.*g1.^2).*tanh(g0)./cosh(g0).^2;
  fr3 = g3./cosh(g0).^2 - 2.*g1.^3 ./ cosh(g0).^4 - ...
         6.*g1.*g2./cosh(g0).^2 .* tanh(g0) + ...
         4.*g1.^3 ./ cosh(g0).^2 .* tanh(g0).^2;
  
  fr = dual3(fr0,fr1,fr2,fr3);
end
 