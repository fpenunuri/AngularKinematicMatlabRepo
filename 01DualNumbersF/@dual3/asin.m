%ASIN   Inverse sine, result in radians.
%   ASIN(X) is the arcsine of the elements of X. Complex
%   results are obtained if ABS(x) > 1.0 for some element.
% The elements of X can be dual3 numbers
function fr = asin(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  
  k = 1 - g0.^2;
  
  fr0 = asin(g0);
  fr1 = g1./sqrt(k);
  fr2 = g2./sqrt(k) + (g0./k.^(3/2)) .* g1.^2;
  fr3 = (3.*g0.^2 .* g1.^3)./(1 - g0.^2).^2.5 + g1.^3./ ...
        (1 - g0.^2).^1.5 + (3.*g0.*g1.*g2)./(1 - g0.^2).^1.5 + ... 
        g3./sqrt(1 - g0.^2);
  
  fr = dual3(fr0,fr1,fr2,fr3);
end