%ATANH  Inverse hyperbolic tangent.
%   ATANH(X) is the inverse hyperbolic tangent of the elements of X.
% The elements of X can be dual4 numbers
function fr = atanh(g)
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  g4 = g.f4;
  
  vd3 = atanh(dual3(g0,g1,g2,g3));  
  vd4 = (6.*(4.*g0.*(1 + g0.^2).*g1.^4 + 2.*(1 + 2.*g0.^2 - 3.*g0.^4).* ...
    g1.^2.*g2 + g0.*(-1 + g0.^2).^2.*g2.^2) + 8.*g0.*(-1 + ...
    g0.^2).^2 .* g1.*g3 - (-1 + g0.^2).^3.*g4)./(-1 + g0.^2).^4;
  
  fr = dual4(vd3.f0,vd3.f1,vd3.f2,vd3.f3,vd4);
end

