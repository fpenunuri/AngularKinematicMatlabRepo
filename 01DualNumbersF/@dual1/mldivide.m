%\   Backslash or left matrix divide.
%   A\B is the matrix division of A into B, which is roughly the
%   same as INV(A)*B 
% dual1 case

function fr = mldivide(A,B)

    A0 = A.f0;
    A1 = A.f1;
    ia0 = A0\eye(size(A0));
    
    fr0 = ia0;
    fr1 = -ia0*A1*ia0;
  
    IA = dual1(fr0,fr1);
    fr=IA*B;  
end



