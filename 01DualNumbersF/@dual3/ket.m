%ketA = ket(A)
%use this function to build a ket vector with the elements of A
%ketA = A(:) only works outside the folder @dual3

%F. Penunuri
%UADY
function fr = ket(x)
    if (isa(x,'dual3'))
        x0 = x.f0;
        x1 = x.f1;
        x2 = x.f2;
        x3 = x.f3;
    
        ket0 = x0(:);
        ket1 = x1(:);
        ket2 = x2(:);
        ket3 = x3(:);

        fr = dual3(ket0,ket1,ket2,ket3);
    else
        fr = x(:);
    end        
end






