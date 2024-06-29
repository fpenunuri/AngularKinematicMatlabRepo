%.' Transpose.
%   X.' is the non-conjugate transpose.
%
%The elements of A can be dual1 numbers
function fr = transpose(A)
  A0 = A.f0;
  A1 = A.f1;

  r0 = builtin('transpose', A0);
  r1 = builtin('transpose', A1);
  
  fr = dual1(r0,r1);
end

