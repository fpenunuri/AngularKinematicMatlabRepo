%use this function to create a new array from A, changing the elements 
%'varargin' for x.
%Example: A=setpart(A,x,i,j,k) --> A(i,j,k) = x
function fr = setpart(A,x,varargin)
  if(isnumeric(A) && isnumeric(x))   
    A(varargin{:}) = x;
    fr = A;
    return
  end
  
  if(isnumeric(A) && isa(x,'dual2'))
    A = dual2(A);
  elseif(isa(A,'dual2') && isnumeric(x))
    x = dual2(x);
  end
  
  if(isa(A,'dual2'))  
    A0 = A.f0; 
    A1 = A.f1; 
    A2 = A.f2; 
    
    x0 = x.f0; 
    x1 = x.f1;
    x2 = x.f2; 
    
    A0(varargin{:}) = x0;
    A1(varargin{:}) = x1;
    A2(varargin{:}) = x2;
    
    fr = dual2(A0,A1,A2);
  else
    S.subs = varargin;
    fr = builtin('subsasgn',A,S,x);
  end
end
