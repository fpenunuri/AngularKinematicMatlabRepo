%ERF Error function.
%   Y = ERF(X) is the error function for each element of X.  X must be
%   real. The error function is defined as:
%
%     erf(x) = 2/sqrt(pi) * integral from 0 to x of exp(-t^2) dt.
%
% The elements of X can be dual2 numbers
function fr = erf(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  
  k = 2.0./sqrt(pi);

  fr0 = erf(g0);
  fr1 = k.*exp(-g0.^2).*g1;
  fr2 = k.*g2.*exp(-g0.^2) - 2.*k.*g0.*exp(-g0.^2).*g1.^2;
  
  fr = dual2(fr0,fr1,fr2);
end
