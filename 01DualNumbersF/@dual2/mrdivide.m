%/   Slash or right matrix divide.
%   A/B is the matrix division of B into A, which is roughly the
%   same as A*INV(B) , except it is computed in a different way.
% dual2 case
function fr = mrdivide(A,B)

  fr = A*inv(B);
end
