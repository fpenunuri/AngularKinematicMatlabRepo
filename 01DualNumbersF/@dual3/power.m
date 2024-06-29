%.^  Array power.
% The elements of X and Y can be dual3 numbers

function fr = power(x,y)
    A=dual3(x);
    B=dual3(y);

    A0 = A.f0; A1 = A.f1; A2 = A.f2; A3 = A.f3;
    B0 = B.f0; B1 = B.f1; B2 = B.f2; B3 = B.f3;

    frd2 = power(dual2(A0,A1,A2),dual2(B0,B1,B2));

    fr3 = A0.^(B0 - 3).*(2.*A1.^3.*B0 - 3.*A0.*A1.*A2.*B0 + A0.^2.*A3.* ...
    B0 - 3.*A0.*A1.^2.*B1 + 3.*A0.^2.*A2.*B1 + 3.*A0.^2.*A1.*B2 + ...
    A0.^3.*B3.*log(A0) + (A1.*B0 + A0.*B1.*log(A0)).^3 - 3.*(A1.* ...
    B0 + A0.*B1.*log(A0)).*(A1.^2.*B0 - A0.*A2.*B0 - 2.*A0.*A1.*B1 - ...
    A0.^2.*B2.*log(A0)));    

    fr3(B0==0 & B1==0 & B2==0 & B3==0) = 0;
    fr3(A0==0 & A1==0 & A2==0 & A3==0 & real(B0)>0) = 0;
    fr3(A1==0 & A2==0 & A3==0 & B1==0 & B2==0 & B3==0) = 0;
    fr3(A0==0 & real(B0)>3) = 0;
    fr3(A0==0 & B1==0 & B2==0 & B0==1) = A3(A0==0 & B1==0 & B2==0 & B0==1);
    fr3(A0==0 & B1==0 & B2==0 & B0==2) = ...
        6.*A1(A0==0 & B1==0 & B2==0 & B0==2).* ...
           A2(A0==0 & B1==0 & B2==0 & B0==2);

    fr3(A0==0 & B1==0 & B2==0 & B0==3) = ...
        6.*A1(A0==0 & B1==0 & B2==0 & B0==3).^3;   

    fr = dual3(frd2.f0,frd2.f1,frd2.f2,fr3);
end

% $$$ function fr = power(x,y)
% $$$     A=dual3(x);
% $$$     B=dual3(y);
% $$$ 
% $$$     A0 = A.f0; A1 = A.f1; A2 = A.f2; A3 = A.f3;
% $$$     B0 = B.f0; B1 = B.f1; B2 = B.f2; B3 = B.f3;
% $$$ 
% $$$     numzero = eps;
% $$$     
% $$$     %this function may need some adjustments
% $$$     A0(A0==0)=numzero; A1(A1==0)=numzero; A2(A2==0)=numzero; A3(A3==0)=numzero;
% $$$     B0(B0==0)=numzero; B1(B1==0)=numzero; B2(B2==0)=numzero; B3(B3==0)=numzero;
% $$$     
% $$$     A = dual3(A0,A1,A2,A3); B = dual3(B0,B1,B2,B3);
% $$$     
% $$$     fr = exp(B.*log(A));
% $$$ end
