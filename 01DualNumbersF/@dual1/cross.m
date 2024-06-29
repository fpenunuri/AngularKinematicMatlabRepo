%cross(X,Y,dim) product
%The elements of X and Y can be dual1 numbers

%F. Penunuri
%UADY
function fr = cross(varargin)
if(nargin==2)
    A = varargin{1};
    B = varargin{2};    
    fr = cross2argsd1(A,B);
elseif(nargin==3)
    A = varargin{1};
    B = varargin{2};    
    dim = varargin{3};        
    fr = cross3argsd1(A,B,dim);        
else
    error('wrong number of arguments in dual1/cross function')          
end
end

function fr = cross3argsd1(A,B,dim)
if(isa(A,'dual1') && isa(B,'dual1'))
    A0=A.f0; A1= A.f1;
    B0=B.f0; B1= B.f1;
    fr0 = cross(A0,B0,dim);
    fr1 = cross(A1,B0,dim) + cross(A0,B1,dim);
    fr  = dual1(fr0,fr1);
elseif(isa(A,'dual1') && ~(isa(B,'dual1')))
    B=dual1(B);
    A0=A.f0; A1= A.f1;
    B0=B.f0; B1= B.f1;    
    fr0 = cross(A0,B0,dim);
    fr1 = cross(A1,B0,dim) + cross(A0,B1,dim);
    fr  = dual1(fr0,fr1);
elseif(~(isa(A,'dual1')) && isa(B,'dual1'))
    A=dual1(A);
    A0=A.f0; A1= A.f1;
    B0=B.f0; B1= B.f1;    
    fr0 = cross(A0,B0,dim);
    fr1 = cross(A1,B0,dim) + cross(A0,B1,dim);
    fr  = dual1(fr0,fr1);
else
    fr = cross(A,B,dim);
end
end

function fr = cross2argsd1(A,B)
if(isa(A,'dual1') && isa(B,'dual1'))
    A0=A.f0; A1= A.f1;
    B0=B.f0; B1= B.f1;
    fr0 = cross(A0,B0);
    fr1 = cross(A1,B0) + cross(A0,B1);
    fr = dual1(fr0,fr1);
elseif(isa(A,'dual1') && ~(isa(B,'dual1')))
    B=dual1(B);
    A0=A.f0; A1= A.f1;
    B0=B.f0; B1= B.f1;    
    fr0 = cross(A0,B0);
    fr1 = cross(A1,B0) + cross(A0,B1);
    fr = dual1(fr0,fr1);
elseif(~(isa(A,'dual1')) && isa(B,'dual1'))
    A=dual1(A);
    A0=A.f0; A1= A.f1;
    B0=B.f0; B1= B.f1;    
    fr0 = cross(A0,B0);
    fr1 = cross(A1,B0) + cross(A0,B1);
    fr = dual1(fr0,fr1);
else
    fr = cross(A,B);
end
end

