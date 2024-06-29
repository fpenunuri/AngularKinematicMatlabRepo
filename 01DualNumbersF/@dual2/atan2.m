%ATAN2  Four quadrant inverse tangent.
%   ATAN2(Y,X) is the four quadrant arctangent of the elements
%   of X and Y. -pi <= ATAN2(Y,X) <= pi.
% The elements of Y and X can be dual2 numbers
function fr = atan2(A,B)
  if (isa(A,'dual2') && isa(B,'dual2'))
    A0 = A.f0;
    A1 = A.f1;
    A2 = A.f2;
    
    B0 = B.f0;
    B1 = B.f1;
    B2 = B.f2;
    
  elseif (isa(A,'dual2') && isnumeric(B))
    A0 = A.f0;
    A1 = A.f1;
    A2 = A.f2;
    
    B0 = B;
    B1 = zeros(size(B));
    B2 = B1;
    
  elseif (isnumeric(A) && isa(B,'dual2'))
    B0 = B.f0;
    B1 = B.f1;
    B2 = B.f2;
    
    A0 = A;
    A1 = zeros(size(A));
    A2 = A1;
    
  end

  r0 = atan2_z(A0,B0);
  r1 = (A1.*B0 - A0.*B1)./(A0.^2 + B0.^2);
  r2 = (B0.^2 .* (A2.*B0 - 2.*A1.*B1) + ...
        A0.^2 .* (A2.*B0 + 2.*A1.*B1) - ...
        A0.^3 .* B2 - A0.*B0.*(2.*A1.^2 - 2.*B1.^2 + B0.*B2))./ ...
        (A0.^2 + B0.^2).^2;

  fr = dual2(r0,r1,r2);
end
