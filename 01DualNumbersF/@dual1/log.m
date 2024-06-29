%LOG    Natural logarithm.
% The elements of X can be dual1 numbers
function fr = log(g)
  g0 = g.f0;
  g1 = g.f1;
  
  fr0 = log(g0);
  fr1 = g1./g0;
  
  fr = dual1(fr0,fr1);
end