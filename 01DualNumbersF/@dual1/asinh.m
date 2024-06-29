%ASINH  Inverse hyperbolic sine.
%   ASINH(X) is the inverse hyperbolic sine of the elements of X.
% The elements of X can be dual1 numbers
function fr = asinh(g) 
  g0 = g.f0;
  g1 = g.f1;
  
  fr0 = asinh(g0);
  fr1 = g1./sqrt(1 + g0.^2);
  
  fr = dual1(fr0,fr1);
end
