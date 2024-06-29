%==  Equal.
%   A == B 
% A and B can be dual1 numbers
function fr = eq(Ax,Bx)
  
  A = dual1(Ax);
  B = dual1(Bx);
  
  A0 = A.f0;
  A1 = A.f1;

  B0 = B.f0;
  B1 = B.f1;

  fr0 = A0==B0;
  fr1 = A1==B1;
  
  fr = fr0 & fr1;  
end
