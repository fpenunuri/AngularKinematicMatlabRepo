%INV    Matrix inverse.
%   INV(X) is the inverse of the square matrix X.
%   A warning message is printed if X is badly scaled or
%   nearly singular.
% The elements of X can be dual4 numbers

function fr = inv(A)
  
  A0 = A.f0;
  A1 = A.f1;
  A2 = A.f2;
  A3 = A.f3;
  A4 = A.f4;
  
  ia0 = A0\eye(size(A0));
  
  fr0 = ia0;
  fr1 = -ia0*A1*ia0;
  fr2 = 2.*ia0*A1*ia0*A1*ia0 - ia0*A2*ia0;
  fr3 = -fr0*A3*fr0 - 3*fr0*A2*fr1 - 3*fr0*A1*fr2;
  fr4 = -fr0*A4*fr0 - 4*fr0*A3*fr1 - 6*fr0*A2*fr2 - 4*fr0*A1*fr3;
  
  fr = dual4(fr0,fr1,fr2,fr3,fr4);
end

