%dot(X,Y,dim) product
%The elements of X and Y can be dual3 numbers
%|A> = |A0> + |A1> eps + |A2> eps2 + |A3>eps3; |B> = ...
%As in the conjg and abs functions, here, the dual parts are (df)* not
%d(f*)
%F. Penunuri
%UADY
function fr = dot(varargin)
if(nargin==2)
    A = varargin{1};
    B = varargin{2};    
    fr = dot2argsd3(A,B);
elseif(nargin==3)
    A = varargin{1};
    B = varargin{2};    
    dim = varargin{3};        
    fr = dot3argsd3(A,B,dim);        
else
    error('wrong number of arguments in dual3/dot function')          
end
end

function fr = dot3argsd3(A,B,dim)
if(isa(A,'dual3') && isa(B,'dual3'))
    A0 = A.f0; A1 = A.f1; A2 = A.f2; A3 = A.f3;
    B0 = B.f0; B1 = B.f1; B2 = B.f2; B3 = B.f3;
    fr0 = dot(A0,B0,dim);
    fr1 = dot(A1,B0,dim) + dot(A0,B1,dim);
    fr2 = dot(A0,B2,dim) + 2.*dot(A1,B1,dim) + dot(A2,B0,dim);
    fr3 = dot(A0,B3,dim) + 3.*dot(A1,B2,dim) + 3.*dot(A2,B1,dim) + ...
      dot(A3,B0,dim);
  
    fr  = dual3(fr0,fr1,fr2,fr3);
elseif(isa(A,'dual3') && ~(isa(B,'dual3')))
    B = dual3(B);
    A0 = A.f0; A1 = A.f1; A2 = A.f2; A3 = A.f3;
    B0 = B.f0; B1 = B.f1; B2 = B.f2; B3 = B.f3;
    fr0 = dot(A0,B0,dim);
    fr1 = dot(A1,B0,dim) + dot(A0,B1,dim);
    fr2 = dot(A0,B2,dim) + 2.*dot(A1,B1,dim) + dot(A2,B0,dim);
    fr3 = dot(A0,B3,dim) + 3.*dot(A1,B2,dim) + 3.*dot(A2,B1,dim) + ...
      dot(A3,B0,dim);
  
    fr  = dual3(fr0,fr1,fr2,fr3);
elseif(~(isa(A,'dual3')) && isa(B,'dual3'))
    A = dual3(A);
    A0 = A.f0; A1 = A.f1; A2 = A.f2; A3 = A.f3;
    B0 = B.f0; B1 = B.f1; B2 = B.f2; B3 = B.f3;
    fr0 = dot(A0,B0,dim);
    fr1 = dot(A1,B0,dim) + dot(A0,B1,dim);
    fr2 = dot(A0,B2,dim) + 2.*dot(A1,B1,dim) + dot(A2,B0,dim);
    fr3 = dot(A0,B3,dim) + 3.*dot(A1,B2,dim) + 3.*dot(A2,B1,dim) + ...
      dot(A3,B0,dim);
  
    fr  = dual3(fr0,fr1,fr2,fr3);
else
    fr = dot(A,B,dim);
end
end

function fr = dot2argsd3(A,B)
if(isa(A,'dual3') && isa(B,'dual3'))
    A0 = A.f0; A1 = A.f1; A2 = A.f2; A3 = A.f3;
    B0 = B.f0; B1 = B.f1; B2 = B.f2; B3 = B.f3;
    fr0 = dot(A0,B0);
    fr1 = dot(A1,B0) + dot(A0,B1);
    fr2 = dot(A0,B2) + 2.*dot(A1,B1) + dot(A2,B0);
    fr3 = dot(A0,B3) + 3.*dot(A1,B2) + 3.*dot(A2,B1) + dot(A3,B0);  
    fr  = dual3(fr0,fr1,fr2,fr3);
elseif(isa(A,'dual3') && ~(isa(B,'dual3')))
    B=dual3(B);
    A0 = A.f0; A1 = A.f1; A2 = A.f2; A3 = A.f3;
    B0 = B.f0; B1 = B.f1; B2 = B.f2; B3 = B.f3;
    fr0 = dot(A0,B0);
    fr1 = dot(A1,B0) + dot(A0,B1);
    fr2 = dot(A0,B2) + 2.*dot(A1,B1) + dot(A2,B0);
    fr3 = dot(A0,B3) + 3.*dot(A1,B2) + 3.*dot(A2,B1) + dot(A3,B0);
    fr  = dual3(fr0,fr1,fr2,fr3);
elseif(~(isa(A,'dual3')) && isa(B,'dual3'))
    A=dual3(A);
    A0 = A.f0; A1 = A.f1; A2 = A.f2; A3 = A.f3;
    B0 = B.f0; B1 = B.f1; B2 = B.f2; B3 = B.f3;
    fr0 = dot(A0,B0);
    fr1 = dot(A1,B0) + dot(A0,B1);
    fr2 = dot(A0,B2) + 2.*dot(A1,B1) + dot(A2,B0);
    fr3 = dot(A0,B3) + 3.*dot(A1,B2) + 3.*dot(A2,B1) + dot(A3,B0);
    fr  = dual3(fr0,fr1,fr2,fr3);
else
    fr = dot(A,B);
end
end
