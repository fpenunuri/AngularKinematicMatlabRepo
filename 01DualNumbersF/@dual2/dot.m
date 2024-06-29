%dot(X,Y,dim) product
%The elements of X and Y can be dual2 numbers

%F. Penunuri
%UADY
function fr = dot(varargin)
if(nargin==2)
    A = varargin{1};
    B = varargin{2};    
    fr = dot2argsd2(A,B);
elseif(nargin==3)
    A = varargin{1};
    B = varargin{2};    
    dim = varargin{3};        
    fr = dot3argsd2(A,B,dim);        
else
    error('wrong number of arguments in dual2/dot function')          
end
end

function fr = dot3argsd2(A,B,dim)
if(isa(A,'dual2') && isa(B,'dual2'))
    A0=A.f0; A1= A.f1; A2=A.f2;
    B0=B.f0; B1= B.f1; B2=B.f2;
    fr0 = dot(A0,B0,dim);
    fr1 = dot(A1,B0,dim) + dot(A0,B1,dim);
    fr2 = dot(A0,B2,dim) + 2.*dot(A1,B1,dim) + dot(A2,B0,dim);
    fr  = dual2(fr0,fr1,fr2);
elseif(isa(A,'dual2') && ~(isa(B,'dual2')))
    B=dual2(B);
    A0=A.f0; A1= A.f1; A2=A.f2;
    B0=B.f0; B1= B.f1; B2=B.f2;    
    fr0 = dot(A0,B0,dim);
    fr1 = dot(A1,B0,dim) + dot(A0,B1,dim);
    fr2 = dot(A0,B2,dim) + 2.*dot(A1,B1,dim) + dot(A2,B0,dim);
    fr  = dual2(fr0,fr1,fr2);
elseif(~(isa(A,'dual2')) && isa(B,'dual2'))
    A=dual2(A);
    A0=A.f0; A1= A.f1; A2=A.f2;
    B0=B.f0; B1= B.f1; B2=B.f2;    
    fr0 = dot(A0,B0,dim);
    fr1 = dot(A1,B0,dim) + dot(A0,B1,dim);
    fr2 = dot(A0,B2,dim) + 2.*dot(A1,B1,dim) + dot(A2,B0,dim);
    fr  = dual2(fr0,fr1,fr2);
else
    fr = dot(A,B,dim);
end
end

function fr = dot2argsd2(A,B)
if(isa(A,'dual2') && isa(B,'dual2'))
    A0=A.f0; A1= A.f1; A2=A.f2;
    B0=B.f0; B1= B.f1; B2=B.f2;
    fr0 = dot(A0,B0);
    fr1 = dot(A1,B0) + dot(A0,B1);
    fr2 = dot(A0,B2) + 2.*dot(A1,B1) + dot(A2,B0);
    fr = dual2(fr0,fr1,fr2);
elseif(isa(A,'dual2') && ~(isa(B,'dual2')))
    B=dual2(B);
    A0=A.f0; A1= A.f1; A2=A.f2;
    B0=B.f0; B1= B.f1; B2=B.f2;    
    fr0 = dot(A0,B0);
    fr1 = dot(A1,B0) + dot(A0,B1);
    fr2 = dot(A0,B2) + 2.*dot(A1,B1) + dot(A2,B0);
    fr = dual2(fr0,fr1,fr2);
elseif(~(isa(A,'dual2')) && isa(B,'dual2'))
    A=dual2(A);
    A0=A.f0; A1= A.f1; A2=A.f2;
    B0=B.f0; B1= B.f1; B2=B.f2;    
    fr0 = dot(A0,B0);
    fr1 = dot(A1,B0) + dot(A0,B1);
    fr2 = dot(A0,B2) + 2.*dot(A1,B1) + dot(A2,B0);
    fr = dual2(fr0,fr1,fr2);
else
    fr = dot(A,B);
end
end

