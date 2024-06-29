%.' Transpose.
%   X.' is the non-conjugate transpose.
%
%The elements of A can be dual3 numbers
function fr = transpose(A)
  A0 = A.f0;
  A1 = A.f1;
  A2 = A.f2;
  A3 = A.f3;

  r0 = builtin('transpose', A0);
  r1 = builtin('transpose', A1);
  r2 = builtin('transpose', A2);
  r3 = builtin('transpose', A3);
  
  fr = dual3(r0,r1,r2,r3);
end

