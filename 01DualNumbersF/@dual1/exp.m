%EXP    Exponential.
% The elements of X can be dual1 numbers
function fr = exp(g) 
  
  g0 = g.f0;
  g1 = g.f1;
  
  fr0 = exp(g0);
  fr1 = exp(g0).*g1;
  
  fr = dual1(fr0,fr1);
end
