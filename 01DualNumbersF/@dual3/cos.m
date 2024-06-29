%COS    Cosine of argument in radians.
% The elements of X can be dual3 numbers
function fr = cos(g)
% fr = sin(g + pi/2);
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  
  cd2 = cos(dual2(g0,g1,g2));  
  c3 = -3.*g1.*g2.*cos(g0) + g1.^3 .* sin(g0) - g3.*sin(g0);
  
  fr = dual3(cd2.f0,cd2.f1,cd2.f2,c3);
end
