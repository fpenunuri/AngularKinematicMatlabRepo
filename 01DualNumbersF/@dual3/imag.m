%IMAG   Complex imaginary part.
%   IMAG(X) is the imaginary part of X.
% dual3 case

function fr = imag(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;  
  
  fr0 = imag(g0);
  fr1 = imag(g1);
  fr2 = imag(g2);
  fr3 = imag(g3);
  
  fr = dual3(fr0,fr1,fr2,fr3);
end
