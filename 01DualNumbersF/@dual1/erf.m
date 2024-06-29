%ERF Error function.
%   Y = ERF(X) is the error function for each element of X.  X must be
%   real. The error function is defined as:
%
%     erf(x) = 2/sqrt(pi) * integral from 0 to x of exp(-t^2) dt.
%
% The elements of X can be dual1 numbers
function fr = erf(g) 
  g0 = g.f0;
  g1 = g.f1;
  
  k = 2.0./sqrt(pi);

  fr0 = erf(g0);
  fr1 = k.*exp(-g0.^2).*g1;
  
  fr = dual1(fr0,fr1);
end
