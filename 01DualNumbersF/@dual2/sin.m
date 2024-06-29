% The elements of X can be dual2 numbers
function local = sin(g)
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;

  sd1 = sin(dual1(g0,g1));  
  sin2 = -sin(g0).*(g1.^2) + cos(g0).*g2;
  
  local = dual2(sd1.f0,sd1.f1,sin2);
end
