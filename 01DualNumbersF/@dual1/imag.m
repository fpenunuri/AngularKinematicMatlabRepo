%IMAG   Complex imaginary part.
%   IMAG(X) is the imaginary part of X.
% dual 1 case

function fr = imag(g) 
  g0 = g.f0;
  g1 = g.f1;
  
  fr0 = imag(g0);
  fr1 = imag(g1);  
  
  fr = dual1(fr0,fr1);
end
