%use this function to create a new array from A, changing the elements 
%'varargin' for x.
%Example: A=setpart(A,x,i,j,k) --> A(i,j,k) = x
% dual4 case 

function fr = setpart(A,x,varargin)
  if(isnumeric(A) && isnumeric(x))   
    A(varargin{:}) = x;
    fr = A;
    return
  end
  
  if(isnumeric(A) && isa(x,'dual4'))
    A = dual4(A);
  elseif(isa(A,'dual4') && isnumeric(x))
    x = dual4(x);
  end
  
  if(isa(A,'dual4'))  
    A0 = A.f0; 
    A1 = A.f1; 
    A2 = A.f2; 
    A3 = A.f3; 
    A4 = A.f4; 
    
    x0 = x.f0; 
    x1 = x.f1;
    x2 = x.f2;
    x3 = x.f3;
    x4 = x.f4;
    
    A0(varargin{:}) = x0;
    A1(varargin{:}) = x1;
    A2(varargin{:}) = x2;
    A3(varargin{:}) = x3;
    A4(varargin{:}) = x4;
    
    fr = dual4(A0,A1,A2,A3,A4);
  else
    S.subs = varargin;
    fr = builtin('subsasgn',A,S,x);
  end
end
