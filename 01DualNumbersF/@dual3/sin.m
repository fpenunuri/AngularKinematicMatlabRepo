% The elements of X can be dual3 numbers
function local = sin(g)
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  
  sd2 = sin(dual2(g0,g1,g2));  
  sin3 = -cos(g0) .* g1.^3 + g3 .* cos(g0) - 3 .* g1 .* g2 .* sin(g0);
  
  local = dual3(sd2.f0,sd2.f1,sd2.f2,sin3);
end
