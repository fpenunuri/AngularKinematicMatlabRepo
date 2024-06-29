%cross(X,Y,dim) product
%The elements of X and Y can be dual2 numbers

%F. Penunuri
%UADY
function fr = cross(varargin)
if(nargin==2)
    A = varargin{1};
    B = varargin{2};    
    fr = cross2argsd2(A,B);
elseif(nargin==3)
    A = varargin{1};
    B = varargin{2};    
    dim = varargin{3};        
    fr = cross3argsd2(A,B,dim);        
else
    error('wrong number of arguments in dual2/cross function')          
end
end

function fr = cross3argsd2(A,B,dim)
if(isa(A,'dual2') && isa(B,'dual2'))
    A0=A.f0; A1= A.f1; A2=A.f2;
    B0=B.f0; B1= B.f1; B2=B.f2;
    fr0 = cross(A0,B0,dim);
    fr1 = cross(A1,B0,dim) + cross(A0,B1,dim);
    fr2 = cross(A0,B2,dim) + 2.*cross(A1,B1,dim) + cross(A2,B0,dim);
    fr  = dual2(fr0,fr1,fr2);
elseif(isa(A,'dual2') && ~(isa(B,'dual2')))
    B=dual2(B);
    A0=A.f0; A1= A.f1; A2=A.f2;
    B0=B.f0; B1= B.f1; B2=B.f2;    
    fr0 = cross(A0,B0,dim);
    fr1 = cross(A1,B0,dim) + cross(A0,B1,dim);
    fr2 = cross(A0,B2,dim) + 2.*cross(A1,B1,dim) + cross(A2,B0,dim);
    fr  = dual2(fr0,fr1,fr2);
elseif(~(isa(A,'dual2')) && isa(B,'dual2'))
    A=dual2(A);
    A0=A.f0; A1= A.f1; A2=A.f2;
    B0=B.f0; B1= B.f1; B2=B.f2;    
    fr0 = cross(A0,B0,dim);
    fr1 = cross(A1,B0,dim) + cross(A0,B1,dim);
    fr2 = cross(A0,B2,dim) + 2.*cross(A1,B1,dim) + cross(A2,B0,dim);
    fr  = dual2(fr0,fr1,fr2);
else
    fr = cross(A,B,dim);
end
end

function fr = cross2argsd2(A,B)
if(isa(A,'dual2') && isa(B,'dual2'))
    A0=A.f0; A1= A.f1; A2=A.f2;
    B0=B.f0; B1= B.f1; B2=B.f2;
    fr0 = cross(A0,B0);
    fr1 = cross(A1,B0) + cross(A0,B1);
    fr2 = cross(A0,B2) + 2.*cross(A1,B1) + cross(A2,B0);
    fr = dual2(fr0,fr1,fr2);
elseif(isa(A,'dual2') && ~(isa(B,'dual2')))
    B=dual2(B);
    A0=A.f0; A1= A.f1; A2=A.f2;
    B0=B.f0; B1= B.f1; B2=B.f2;    
    fr0 = cross(A0,B0);
    fr1 = cross(A1,B0) + cross(A0,B1);
    fr2 = cross(A0,B2) + 2.*cross(A1,B1) + cross(A2,B0);
    fr = dual2(fr0,fr1,fr2);
elseif(~(isa(A,'dual2')) && isa(B,'dual2'))
    A=dual2(A);
    A0=A.f0; A1= A.f1; A2=A.f2;
    B0=B.f0; B1= B.f1; B2=B.f2;    
    fr0 = cross(A0,B0);
    fr1 = cross(A1,B0) + cross(A0,B1);
    fr2 = cross(A0,B2) + 2.*cross(A1,B1) + cross(A2,B0);
    fr = dual2(fr0,fr1,fr2);
else
    fr = cross(A,B);
end
end

