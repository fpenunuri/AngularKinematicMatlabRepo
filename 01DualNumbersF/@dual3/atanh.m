%ATANH  Inverse hyperbolic tangent.
%   ATANH(X) is the inverse hyperbolic tangent of the elements of X.
% The elements of X can be dual3 numbers
function fr = atanh(g)
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  
  %Matlab (2020a) doesn't give the principal value for atanh
  fr0 = 1/2*(log(1 + g0) - log(1-g0)); %atanh(g0); 
  fr1 = g1./(1 - g0.^2);
  fr2 = 2.*g0.*g1.^2 ./ (1 - g0.^2).^2 + g2./(1 - g0.^2);
  fr3 = (8.*g0.^2.*g1.^3)./(1 - g0.^2).^3 + (2.*g1.^3)./ ...
        (1 - g0.^2).^2 + (6.*g0.*g1.*g2)./(1 - g0.^2).^2 + ...
        g3./(1 - g0.^2);

  fr = dual3(fr0,fr1,fr2,fr3);
end

