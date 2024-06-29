%ATAN2  Four quadrant inverse tangent.
%   ATAN2(Y,X) is the four quadrant arctangent of the elements
%   of X and Y. -pi <= ATAN2(Y,X) <= pi.
% The elements of Y and X can be dual1 numbers

function fr = atan2(A,B)
  if (isa(A,'dual1') && isa(B,'dual1'))
    A0 = A.f0;
    A1 = A.f1;
    
    B0 = B.f0;
    B1 = B.f1;
    
  elseif (isa(A,'dual1') && isnumeric(B))
    A0 = A.f0;
    A1 = A.f1;
    
    B0 = B;
    B1 = zeros(size(B));
    
  elseif (isnumeric(A) && isa(B,'dual1'))
    B0 = B.f0;
    B1 = B.f1;
    
    A0 = A;
    A1 = zeros(size(A));
    
  end

  r0 = atan2_z(A0,B0);
  r1 = (A1.*B0 - A0.*B1)./(A0.^2 + B0.^2);

  fr = dual1(r0,r1);
end


  
