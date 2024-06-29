%ACOSH  Inverse hyperbolic cosine.
%   ACOSH(X) is the inverse hyperbolic cosine of the elements of X.
%The elements of X can be dual1 numbers
function fr = acosh(g) 
  g0 = g.f0;
  g1 = g.f1;
  
  g0m   = g0 - 1;
  g0p   = g0 + 1;
  sg0m  = sqrt(g0m); 
  sg0p  = sqrt(g0p);

  fr0 = acosh(g0);
  fr1 = g1./(sg0m.*sg0p);
  
  fr = dual1(fr0,fr1);
end

