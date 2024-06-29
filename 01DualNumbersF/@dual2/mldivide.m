%\   Backslash or left matrix divide.
%   A\B is the matrix division of A into B, which is roughly the
%   same as INV(A)*B 
% dual2 case
function fr = mldivide(A,B)

  fr = inv(A)*B;
end
