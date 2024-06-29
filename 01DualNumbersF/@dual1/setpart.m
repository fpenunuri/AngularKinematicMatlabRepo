%use this function to create an array from A, changing the elements 
%'varargin' for x.
%Example: A=setpart(A,x,i,j,k) --> A(i,j,k) = x
function fr = setpart(A,x,varargin)
  if(isnumeric(A) && isnumeric(x))   
    A(varargin{:}) = x;
    fr = A;
    return
  end
  
  if(isnumeric(A) && isa(x,'dual1'))
    A = dual1(A);
  elseif(isa(A,'dual1') && isnumeric(x))
    x = dual1(x);
  end
  
  if(isa(A,'dual1'))  
    A0 = A.f0; 
    A1 = A.f1; 
    
    x0 = x.f0; 
    x1 = x.f1;
    
    A0(varargin{:}) = x0;
    A1(varargin{:}) = x1;
    
    fr = dual1(A0,A1);
  else
    S.subs = varargin;
    fr = builtin('subsasgn',A,S,x);
  end
end
