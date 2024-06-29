%use this function to create a new array from A, changing the elements 
%'varargin' for x.
%Example: A=setpart(A,x,i,j,k) --> A(i,j,k) = x
function fr = setpart(A,x,varargin)
  if(isnumeric(A) && isnumeric(x))   
    A(varargin{:}) = x;
    fr = A;
    return
  end
  
  if(isnumeric(A) && isa(x,'dual3'))
    A = dual3(A);
  elseif(isa(A,'dual3') && isnumeric(x))
    x = dual3(x);
  end
  
  if(isa(A,'dual3'))  
    A0 = A.f0; 
    A1 = A.f1; 
    A2 = A.f2; 
    A3 = A.f3; 
    
    x0 = x.f0; 
    x1 = x.f1;
    x2 = x.f2;
    x3 = x.f3;
    
    A0(varargin{:}) = x0;
    A1(varargin{:}) = x1;
    A2(varargin{:}) = x2;
    A3(varargin{:}) = x3;
    
    fr = dual3(A0,A1,A2,A3);
  else
    S.subs = varargin;
    fr = builtin('subsasgn',A,S,x);
  end
end
