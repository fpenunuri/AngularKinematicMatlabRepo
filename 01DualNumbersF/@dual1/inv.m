%INV    Matrix inverse.
%   INV(X) is the inverse of the square matrix X.
%   A warning message is printed if X is badly scaled or
%   nearly singular.
% The elements of X can be dual1 numbers
% 
function fr = inv(A)
  
  A0 = A.f0;
  A1 = A.f1;
  
  %    ia0 = inv(A0);
  ia0 = A0\eye(size(A0));
  
  fr0 = ia0;
  fr1 = -ia0*A1*ia0;
  
  fr = dual1(fr0,fr1);
end


