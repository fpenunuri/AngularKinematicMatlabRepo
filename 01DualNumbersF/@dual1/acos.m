%ACOS   Inverse cosine, result in radians.
%   ACOS(X) is the arccosine of the elements of X. Complex
%   results are obtained if ABS(x) > 1.0 for some element.
% The elements of X can be dual1 numbers
function fr = acos(g)
  g0 = g.f0;
  g1 = g.f1;
  
  kd = 1 - g0.^2;
  skd = sqrt(kd);
  
  fr0 = acos(g0);
  fr1 = -g1./skd;
  
  fr = dual1(fr0,fr1);
end
