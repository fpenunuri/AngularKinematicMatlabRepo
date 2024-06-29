%EXP    Exponential.
% The elements of X can be dual2 numbers
function fr = exp(g) 
  
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  
  fr0 = exp(g0);
  fr1 = exp(g0).*g1;
  fr2 = exp(g0).*(g1.^2) + exp(g0).*g2;
  
  fr = dual2(fr0,fr1,fr2);
end
