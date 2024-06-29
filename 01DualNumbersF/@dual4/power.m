%.^  Array power.
% The elements of X and Y can be dual4 numbers

function fr = power(x,y)
    A=dual4(x);
    B=dual4(y);

    A0 = A.f0; A1 = A.f1; A2 = A.f2; A3 = A.f3; A4 = A.f4;
    B0 = B.f0; B1 = B.f1; B2 = B.f2; B3 = B.f3; B4 = B.f4;

    frd3 = power(dual3(A0,A1,A2,A3),dual3(B0,B1,B2,B3));

    fr4 = A0.^(B0 - 4).*(-6.*A1.^4.*B0 + 12.*A0.*A1.^2.*A2.*B0 - 3.*A0.^2.* ...
    A2.^2.*B0 - 4.*A0.^2.*A1.*A3.*B0 + A0.^3.*A4.*B0 + ...
    8.*A0.*A1.^3.*B1 - 12.*A0.^2.*A1.*A2.*B1 + 4.*A0.^3.*A3.*B1 - ...
    6.*A0.^2.*A1.^2.*B2 + 6.* A0.^3.*A2.*B2 + 4.*A0.^3.*A1.*B3 + ...
    A0.^4.*B4.*log(A0) + (A1.*B0 + A0.*B1.*log(A0)).^4 - 6.*(A1.* ...
    B0 + A0.*B1.*log(A0)).^2.*(A1.^2.*B0 - A0.*A2.*B0 - 2.*A0.*A1 .* ...
    B1 - A0.^2.*B2.*log(A0)) + 3.*(-(A1.^2.*B0)+ A0.*A2.*B0 + ...
    2.*A0.*A1.*B1 + A0.^2.*B2.*log(A0)).^2 + 4.*(A1.*B0 + A0.*B1.* ...
    log(A0)).*(2.*A1.^3.*B0 - 3.*A0.*A1.^2.*B1 + A0.^2.*(A3.*B0 + ...
    3.*A2.*B1) + 3.*A0.*A1.*(-(A2.*B0) + A0.*B2) + A0.^3.*B3.*log(A0)));
    

    fr4(B0==0 & B1==0 & B2==0 & B3==0 & B4==0) = 0;
    fr4(A0==0 & A1==0 & A2==0 & A3==0 & A4==0 & real(B0)>0) = 0;
    fr4(A1==0 & A2==0 & A3==0 & A4==0 & B1==0 & B2==0 & B3==0 &B4==0) = 0;
    fr4(A0==0 & real(B0)>4) = 0;
    fr4(A0==0 & B1==0 & B2==0 & B3==0 & B0==1) = ...
           A4(A0==0 & B1==0 & B2==0 & B3==0 & B0==1);
    
    fr4(A0==0 & B1==0 & B2==0 & B3==0 & B0==2) = ...
           6.*A2(A0==0 & B1==0 & B2==0 & B3==0 & B0==2).^2 + ...
           8.*A1(A0==0 & B1==0 & B2==0 & B3==0 & B0==2).*...
              A3(A0==0 & B1==0 & B2==0 & B3==0 & B0==2);
           

    fr4(A0==0 & B1==0 & B2==0 & B3==0 & B0==3) = ...
           36.*A1(A0==0 & B1==0 & B2==0 & B3==0 & B0==3).^2 .* ...
           A2(A0==0 & B1==0 & B2==0 & B3==0 & B0==3);

    fr4(A0==0 & B1==0 & B2==0 & B3==0 & B0==4) = ...
           24.*A1(A0==0 & B1==0 & B2==0 & B3==0 & B0==4).^4;

    fr = dual4(frd3.f0,frd3.f1,frd3.f2,frd3.f3,fr4);
end

% $$$ 
% $$$ function fr = power(x,y)
% $$$     A=dual4(x);
% $$$     B=dual4(y);
% $$$ 
% $$$     A0 = A.f0; A1 = A.f1; A2 = A.f2; A3 = A.f3; A4 = A.f4;
% $$$     B0 = B.f0; B1 = B.f1; B2 = B.f2; B3 = B.f3; B4 = B.f4;
% $$$ 
% $$$     numzero = eps;
% $$$     
% $$$     %this function may need some adjustments
% $$$     A0(A0==0)=numzero; A1(A1==0)=numzero;
% $$$     A2(A2==0)=numzero; A3(A3==0)=numzero;
% $$$     A4(A4==0)=numzero;
% $$$ 
% $$$     B0(B0==0)=numzero; B1(B1==0)=numzero;
% $$$     B2(B2==0)=numzero; B3(B3==0)=numzero;
% $$$     B4(B4==0)=numzero;
% $$$     
% $$$     A = dual4(A0,A1,A2,A3,A4); B = dual4(B0,B1,B2,B3,B4);
% $$$     
% $$$     fr = exp(B.*log(A));
% $$$ end
