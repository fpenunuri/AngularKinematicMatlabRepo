%==  Equal.
%   A == B 
% A and B can be dual4 numbers
function fr = eq(Ax,Bx)
  
  A = dual4(Ax);
  B = dual4(Bx);
  
  A0 = A.f0;
  A1 = A.f1;
  A2 = A.f2;
  A3 = A.f3;  
  A4 = A.f4;  

  B0 = B.f0;
  B1 = B.f1;
  B2 = B.f2;
  B3 = B.f3;
  B4 = B.f4;

  fr0 = A0==B0;
  fr1 = A1==B1;
  fr2 = A2==B2;
  fr3 = A3==B3;
  fr4 = A4==B4;
  
  fr = fr0 & fr1 & fr2 & fr3 & fr4;  
end
