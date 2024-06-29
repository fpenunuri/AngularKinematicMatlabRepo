%ketA = ket(A)
%use this function to build a ket vector with the elements of A
%ketA = A(:) only works outside the folder @dual1

%F. Penunuri
%UADY
function fr = ket(x)
    if (isa(x,'dual1'))
        x0 = x.f0;
        x1 = x.f1;
    
        ket0 = x0(:);
        ket1 = x1(:);

        fr = dual1(ket0,ket1);
    else
        fr = x(:);
    end        
end






