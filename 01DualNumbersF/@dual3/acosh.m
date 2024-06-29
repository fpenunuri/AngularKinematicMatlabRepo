%ACOSH  Inverse hyperbolic cosine.
%   ACOSH(X) is the inverse hyperbolic cosine of the elements of X.
%The elements of X can be dual3 numbers
function fr = acosh(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;  
  
  g0m   = g0 - 1;
  g0p   = g0 + 1;
  sg0m  = sqrt(g0m); 
  sg0p  = sqrt(g0p);

  fr0 = acosh(g0);
  fr1 = g1./(sg0m.*sg0p);
  fr2 = g2./(sg0m.*sg0p) - (g0.*g1.^2)./(g0m.^(3/2).*g0p.^(3/2));
  fr3 = (3.*g1.^3)./(4.*sg0m.*g0p.^2.5) + g1.^3./(2.*g0m.^1.5.* ...
      g0p.^1.5) + (3.*g1.^3)./(4.*g0m.^2.5 .* sg0p) - ...
      (3.*g1.*g2)./(2 .* sg0m.*g0p.^1.5) - (3.*g1.*g2)./ ...
      (2.*g0m.^1.5.*sg0p) + g3./(sg0m.*sg0p);  
  
  fr = dual3(fr0,fr1,fr2,fr3);
end

