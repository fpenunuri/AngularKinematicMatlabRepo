%use this function to select elements, from a dual3 number 
%with arrays components
%Example: takepart(A,1,2) = A(1,2)
function fr = takepart(A,varargin)
  if(isa(A,'dual3'))        
    A0 = A.f0;
    A1 = A.f1; 
    A2 = A.f2; 
    A3 = A.f3; 
    
    r0 = A0(varargin{:});
    r1 = A1(varargin{:});
    r2 = A2(varargin{:});
    r3 = A3(varargin{:});

    fr = dual3(r0,r1,r2,r3);
  else
    S.subs = varargin;
    fr = builtin('subsref',A,S);
  end
end

