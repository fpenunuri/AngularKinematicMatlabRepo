%==  Equal.
%   A == B 
% A and B can be dual3 numbers
function fr = eq(Ax,Bx)
  
  A = dual3(Ax);
  B = dual3(Bx);
  
  A0 = A.f0;
  A1 = A.f1;
  A2 = A.f2;
  A3 = A.f3;  

  B0 = B.f0;
  B1 = B.f1;
  B2 = B.f2;
  B3 = B.f3;

  fr0 = A0==B0;
  fr1 = A1==B1;
  fr2 = A2==B2;
  fr3 = A3==B3;
  
  fr = fr0 & fr1 & fr2 & fr3;  
end
