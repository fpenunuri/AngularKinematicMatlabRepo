%ERF Error function.
%   Y = ERF(X) is the error function for each element of X.  X must be
%   real. The error function is defined as:
%
%     erf(x) = 2/sqrt(pi) * integral from 0 to x of exp(-t^2) dt.
%
% The elements of X can be dual4 numbers

function fr = erf(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  g4 = g.f4;
  
  vd3 = erf(dual3(g0,g1,g2,g3));  
  vd4 = 2./sqrt(pi).*(-8.*g0.^3.*g1.^4 - 12.*g1.^2.*g2 + 24.*g0.^2 .* ...
    g1.^2 .* g2 + 2.*g0.*(6.*g1.^4 - 3.*g2.^2 - 4.*g1.*g3) + ...
    g4)./exp(g0.^2);         
  
  fr = dual4(vd3.f0,vd3.f1,vd3.f2,vd3.f3,vd4);
end
