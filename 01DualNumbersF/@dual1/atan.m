%ATAN   Inverse tangent, result in radians.
%   ATAN(X) is the arctangent of the elements of X.
% The elements of X can be dual1 numbers
function fr = atan(g) 
  g0 = g.f0;
  g1 = g.f1;
  
  fr0 = atan(g0);
  fr1 = g1./(1 + g0.^2);
  
  fr = dual1(fr0,fr1);
end

