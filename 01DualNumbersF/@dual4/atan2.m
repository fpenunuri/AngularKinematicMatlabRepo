%ATAN2  Four quadrant inverse tangent.
%   ATAN2(Y,X) is the four quadrant arctangent of the elements
%   of X and Y. -pi <= ATAN2(Y,X) <= pi.
% The elements of Y and X can be dual4 numbers

function fr = atan2(A,B)
  if (isa(A,'dual4') && isa(B,'dual4'))
    A0 = A.f0;
    A1 = A.f1;
    A2 = A.f2;
    A3 = A.f3;
    A4 = A.f4;
    
    B0 = B.f0;
    B1 = B.f1;
    B2 = B.f2;
    B3 = B.f3;
    B4 = B.f4;

  elseif (isa(A,'dual4') && isnumeric(B))
    A0 = A.f0;
    A1 = A.f1;
    A2 = A.f2;
    A3 = A.f3;
    A4 = A.f4;
    
    B0 = B;
    B1 = zeros(size(B));
    B2 = B1;
    B3 = B1;
    B4 = B1;
    
  elseif (isnumeric(A) && isa(B,'dual4'))
    B0 = B.f0;
    B1 = B.f1;
    B2 = B.f2;
    B3 = B.f3;
    B4 = B.f4;
    
    A0 = A;
    A1 = zeros(size(A));
    A2 = A1;
    A3 = A1;
    A4 = A1;
  end

  r0 = atan2_z(A0,B0);
  r1 = (A1.*B0 - A0.*B1)./(A0.^2 + B0.^2);
  r2 = (B0.^2 .* (A2.*B0 - 2.*A1.*B1) + ...
        A0.^2 .* (A2.*B0 + 2.*A1.*B1) - ...
        A0.^3 .* B2 - A0.*B0.*(2.*A1.^2 - 2.*B1.^2 + B0.*B2))./ ...
        (A0.^2 + B0.^2).^2;

  rr2 = (A0.^2 + B0.^2);

  r3 = (8.*(A1.*B0 - A0.*B1).*(A0.*A1 + B0.*B1).^2 + 2.*(-(A1.^3.*B0) + ...
         A0.*A1.^2.*B1 + A1.*(-3.*A0.*A2.*B0 + 2.*A0.^2 .* B2 - ...
         B0.*(B1.^2 + B0.*B2)) + B1.*(A0.^2 .* A2 - 2.*A2.*B0.^2 + ...
         A0.*(B1.^2 + 3.*B0.*B2))) .* rr2 + ...
        (A3.*B0 + A2.*B1 - A1.*B2 - A0.*B3).*rr2.^2)./rr2.^3;

  r4 = (48.*(-(A1.*B0) + A0.*B1).*(A0.*A1 + B0.*B1).^3 - 24.*(A0.*A1 + ...
         B0.*B1).*(-(A1.^3.*B0) - 2.*A0.*A1.*A2.*B0 + A0.*A1.^2 .* ...
         B1 + A0.^2.*A2.*B1 - A2.*B0.^2.*B1 + A0.^2.*A1.*B2 - ...
         A1.*B0.*(B1.^2 + B0.*B2) + A0.*B1.*(B1.^2 + 2.*B0.*B2)).* ...
         rr2 + 2.*(A1.^2.*(-6.*A2.*B0 + 6.*A0.*B2) + A0.^2.* ...
         (A3.*B1 + 3.*A2.*B2) - 3.*B0.*(A3.*B0.*B1 + 2.*A2.*B1.^2 + ...
         A2.*B0.*B2) + A1.*(-4.*A0.*A3.*B0 + 3.*A0.^2.*B3 - ...
         B0.^2.*B3) + A0.*(-3.*A2.^2.*B0 + 6.*B1.^2.*B2 + ...
         3.*B0.*B2.^2 + 4.*B0.*B1.*B3)).*rr2.^2 + (A4.*B0 + ...
         2.*A3.*B1 - 2.*A1.*B3 - A0.*B4).*rr2.^3)./rr2.^4;

  fr = dual4(r0,r1,r2,r3,r4);
end
