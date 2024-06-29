%==  Equal.
%   A == B 
% A and B can be dual2 numbers
function fr = eq(Ax,Bx)
  
  A = dual2(Ax);
  B = dual2(Bx);
  
  A0 = A.f0;
  A1 = A.f1;
  A2 = A.f2;

  B0 = B.f0;
  B1 = B.f1;
  B2 = B.f2;

  fr0 = A0==B0;
  fr1 = A1==B1;
  fr2 = A2==B2;
  
  fr = fr0 & fr1 & fr2;
end
