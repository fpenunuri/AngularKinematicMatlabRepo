%ACOS   Inverse cosine, result in radians.
%   ACOS(X) is the arccosine of the elements of X. Complex
%   results are obtained if ABS(x) > 1.0 for some element.
% The elements of X can be dual3 numbers
function fr = acos(g)
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  
  kd = 1 - g0.^2;
  skd = sqrt(kd);
  kd32 = kd.^(3./2);
  kd52 = kd.^(5./2);
  
  fr0 = acos(g0);
  fr1 = -g1./skd;
  fr2 = -g2./skd - (g0./kd32).*g1.^2;
  fr3 = -(g1.^3./kd32) - (3.*g0.*g1.*g2)./kd32 - ...
        (3.*g0.^2.*g1.^3)./kd52 - g3./skd;
  
  fr = dual3(fr0,fr1,fr2,fr3);
end
