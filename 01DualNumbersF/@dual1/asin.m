%ASIN   Inverse sine, result in radians.
%   ASIN(X) is the arcsine of the elements of X. Complex
%   results are obtained if ABS(x) > 1.0 for some element.
% The elements of X can be dual1 numbers
function fr = asin(g) 
  g0 = g.f0;
  g1 = g.f1;
  
  k = 1 - g0.^2;
  
  fr0 = asin(g0);
  fr1 = g1./sqrt(k);
  
  fr = dual1(fr0,fr1);
end