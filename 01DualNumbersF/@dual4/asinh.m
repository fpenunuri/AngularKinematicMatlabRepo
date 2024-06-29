%ASINH  Inverse hyperbolic sine.
%   ASINH(X) is the inverse hyperbolic sine of the elements of X.
% The elements of X can be dual4 numbers
function fr = asinh(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  g4 = g.f4;
  
  vd3 = asinh(dual3(g0,g1,g2,g3));  

  vd4 = (-6.*g1.^2 .* g2 + g0.*(9.*g1.^4 - 3.*g2.^2 - 4.*g1.*g3) - ...
    g0.^5 .* (3.*g2.^2 + 4.*g1.*g3) - 2.*g0.^3 .* (3.*(g1.^4 + ...
    g2.^2) + 4.*g1.*g3) + g4 + g0.^6 .* g4 + 3.*g0.^2 .* ...
    (2.*g1.^2 .* g2 + g4) + 3.*g0.^4 .* (4.*g1.^2 .* g2 + g4))./ ...
    (1 + g0.^2).^3.5;
    
  fr = dual4(vd3.f0,vd3.f1,vd3.f2,vd3.f3,vd4);
end
