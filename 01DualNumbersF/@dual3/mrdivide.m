%/   Slash or right matrix divide.
%   A/B is the matrix division of B into A, which is roughly the
%   same as A*INV(B) 
% dual3 case
function fr = mrdivide(A,B)
  fr = A*inv(B);
end
