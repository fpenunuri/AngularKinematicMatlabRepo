%ATANH  Inverse hyperbolic tangent.
%   ATANH(X) is the inverse hyperbolic tangent of the elements of X.
% The elements of X can be dual2 numbers
function fr = atanh(g)
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  
  %Matlab (2020a) doesn't give the principal value for atanh
  fr0 = 1/2*(log(1 + g0) - log(1-g0)); %atanh(g0); 
  fr1 = g1./(1 - g0.^2);
  fr2 = 2.*g0.*g1.^2 ./ (1 - g0.^2).^2 + g2./(1 - g0.^2);
  
  fr = dual2(fr0,fr1,fr2);
end

