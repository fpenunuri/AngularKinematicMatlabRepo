%ATAN   Inverse tangent, result in radians.
%   ATAN(X) is the arctangent of the elements of X.
% The elements of X can be dual3 numbers
function fr = atan(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  
  fr0 = atan(g0);
  fr1 = g1./(1 + g0.^2);
  fr2 = g2./(1 + g0.^2) - (2.*g0./(g0.^2 + 1).^2).*g1.^2;
  fr3 = (8.*g0.^2 .* g1.^3)./(1 + g0.^2).^3 - (2.*g1.^3)./ ...
    (1 + g0.^2).^2 - (6.*g0.*g1.*g2)./(1 + g0.^2).^2 + g3./(1 + g0.^2);
  
  fr = dual3(fr0,fr1,fr2,fr3);
end

