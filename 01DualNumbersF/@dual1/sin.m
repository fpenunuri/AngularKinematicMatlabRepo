% The elements of X can be dual1 numbers
function local = sin(g) 
  g0 = g.f0;
  g1 = g.f1;
  
  sin0 = sin(g0);
  sin1 = cos(g0).*g1;
  
  local = dual1(sin0,sin1);
end
