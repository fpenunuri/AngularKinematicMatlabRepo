%IMAG   Complex imaginary part.
%   IMAG(X) is the imaginary part of X.
% dual4 case 

function fr = imag(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;  
  g4 = g.f4;  
  
  fr0 = imag(g0);
  fr1 = imag(g1);
  fr2 = imag(g2);
  fr3 = imag(g3);
  fr4 = imag(g4);
  
  fr = dual4(fr0,fr1,fr2,fr3,fr4);
end
