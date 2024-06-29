%TANH   Hyperbolic tangent.
% The elements of X can be dual1 numbers
function fr = tanh(g) 
  g0 = g.f0;
  g1 = g.f1;
  
  fr0 = tanh(g0);
  fr1 = g1./cosh(g0).^2;
  
  fr = dual1(fr0,fr1);
end
