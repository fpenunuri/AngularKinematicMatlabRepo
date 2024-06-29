% The elements of X can be dual4 numbers
function local = sin(g)
  g0 = g.f0; g1 = g.f1; g2 = g.f2; g3 = g.f3; g4 = g.f4;
  
  sd3 = sin(dual3(g0,g1,g2,g3));  
  sin4 = -6.*g1.^2 .* g2.*cos(g0) + g4.*cos(g0) + g1.^4 .* sin(g0) - ...
         3.*g2.^2 .* sin(g0) - 4.*g1.*g3.*sin(g0);

  local = dual4(sd3.f0,sd3.f1,sd3.f2,sd3.f3,sin4);
end
