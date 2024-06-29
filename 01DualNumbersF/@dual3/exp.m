%EXP    Exponential.
% The elements of X can be dual3 numbers
function fr = exp(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  
  ed2 = exp(dual2(g0,g1,g2));  
  exp3 =  exp(g0) .* g1.^3 + 3 .* exp(g0) .* g1 .* g2 + exp(g0) .* g3;
  
  fr = dual3(ed2.f0,ed2.f1,ed2.f2,exp3);
end
