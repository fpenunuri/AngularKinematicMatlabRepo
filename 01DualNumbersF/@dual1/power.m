%.^  Array power.
% The elements of X and Y can be dual1 numbers

function fr = power(x,y)
    
    A=dual1(x);
    B=dual1(y);

    A0 = A.f0; A1 = A.f1; B0 = B.f0; B1 = B.f1;

    fr0= A0.^B0;
    fr1 =  A0.^(B0 - 1).*A1.*B0 + A0.^B0.*B1.*log(A0);

    fr1(B0==0 & B1==0) = 0;
    fr1(A0==0 & A1==0 & real(B0)>0) = 0;
    fr1(A1==0 & B1==0) = 0;
    fr1(A0==0 & real(B0)>1) = 0;
    fr1(A0==0 & B0==1) = A1(A0==0 & B0==1);
    
    fr = dual1(fr0,fr1);
end


% $$$ 
% $$$ function fr = power(x,y)
% $$$     A=dual1(x);
% $$$     B=dual1(y);
% $$$ 
% $$$     A0 = A.f0; A1 = A.f1; B0 = B.f0; B1 = B.f1;
% $$$ 
% $$$     numzero = eps;
% $$$     
% $$$     %this function may need some adjustments
% $$$     A0(A0==0)=numzero; A1(A1==0)=numzero; B0(B0==0)=numzero; B1(B1==0)=numzero;
% $$$     
% $$$     A = dual1(A0,A1); B = dual1(B0,B1);
% $$$     
% $$$     fr = exp(B.*log(A));
% $$$ end
