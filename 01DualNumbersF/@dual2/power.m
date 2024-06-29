%.^  Array power.
% The elements of X and Y can be dual2 numbers


function fr = power(x,y)
    A=dual2(x);
    B=dual2(y);

    A0 = A.f0; A1 = A.f1; A2 = A.f2; B0 = B.f0; B1 = B.f1; B2 = B.f2;

    frd1 = power(dual1(A0,A1),dual1(B0,B1));

    fr2 = A0.^B0.*((A1.*B0)./A0 + B1.*log(A0)).^2 + ...
            A0.^B0.*(-((A1.^2.*B0)./A0.^2) + (A2.*B0)./A0 + ... 
            (2.*A1.*B1)./A0 + B2.*log(A0));

    fr2(B0==0 & B1==0 & B2==0) = 0;
    fr2(A0==0 & A1==0 & A2==0 & real(B0)>0) = 0;
    fr2(A1==0 & A2==0 & B1==0 & B2==0) = 0;
    fr2(A0==0 & real(B0)>2) = 0;
    fr2(A0==0 & B1==0 & B0==1) = A2(A0==0 & B1==0 & B0 ==1);
    fr2(A0==0 & B1==0 & B0==2) = 2.*A1(A0==0 & B1==0 & B0 ==2).^2;

    fr = dual2(frd1.f0,frd1.f1,fr2);
end


% $$$ 
% $$$ function fr = power(x,y)
% $$$     A=dual2(x);
% $$$     B=dual2(y);
% $$$ 
% $$$     A0 = A.f0; A1 = A.f1; A2 = A.f2; 
% $$$     B0 = B.f0; B1 = B.f1; B2 = B.f2; 
% $$$ 
% $$$     numzero = eps;
% $$$     
% $$$     %this function may need some adjustments
% $$$     A0(A0==0)=numzero; A1(A1==0)=numzero; A2(A2==0)=numzero;
% $$$     B0(B0==0)=numzero; B1(B1==0)=numzero; B2(B2==0)=numzero;
% $$$     
% $$$     A = dual2(A0,A1,A2); B = dual2(B0,B1,B2);
% $$$     
% $$$     fr = exp(B.*log(A));
% $$$ end
% $$$ 
