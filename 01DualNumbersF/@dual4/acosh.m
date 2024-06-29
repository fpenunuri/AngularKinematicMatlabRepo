%ACOSH  Inverse hyperbolic cosine.
%   ACOSH(X) is the inverse hyperbolic cosine of the elements of X.
%The elements of X can be dual4 numbers

function fr = acosh(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;  
  g4 = g.f4;  
  
  vd3 = acosh(dual3(g0,g1,g2,g3));  

  vd4 = (-9.*g0.*g1.^4 - 6.*g0.^3 .* g1.^4 - 6.*g1.^2 .* g2 - ...
      6.*g0.^2 .* g1.^2 .* g2 + 12.*g0.^4 .* g1.^2 .* g2 - ...
      3.*g0.*g2.^2 + 6.*g0.^3 .* g2.^2 - 3.*g0.^5 .* g2.^2 - ...
      4.*g0.*g1.*g3 + 8.*g0.^3 .* g1.*g3 - 4.*g0.^5 .* g1.*g3 + ...
      (-1 + g0.^2).^3 .* g4)./((-1 + g0).^3.5 .* (1 + g0).^3.5);
  
  fr = dual4(vd3.f0,vd3.f1,vd3.f2,vd3.f3,vd4);
end

