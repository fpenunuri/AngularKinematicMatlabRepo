%ASINH  Inverse hyperbolic sine.
%   ASINH(X) is the inverse hyperbolic sine of the elements of X.
% The elements of X can be dual2 numbers
function fr = asinh(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  
  fr0 = asinh(g0);
  fr1 = g1./sqrt(1 + g0.^2);
  fr2 = g2./sqrt(1 + g0.^2) - (g0.*g1.^2)./(1 + g0.^2).^(3/2);
  
  fr = dual2(fr0,fr1,fr2);
end
