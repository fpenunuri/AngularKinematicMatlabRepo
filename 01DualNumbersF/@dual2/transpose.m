%.' Transpose.
%   X.' is the non-conjugate transpose.
%
%The elements of A can be dual2 numbers
function fr = transpose(A)
  A0 = A.f0;
  A1 = A.f1;
  A2 = A.f2;

  r0 = builtin('transpose', A0);
  r1 = builtin('transpose', A1);
  r2 = builtin('transpose', A2);
  
  fr = dual2(r0,r1,r2);
end

