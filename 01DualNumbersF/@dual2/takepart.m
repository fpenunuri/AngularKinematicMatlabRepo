%use this function to select elements, from a dual2 number 
%with arrays components
%Example: takepart(A,1,2) = A(1,2)
function fr = takepart(A,varargin)
  if(isa(A,'dual2'))        
    A0 = A.f0;
    A1 = A.f1; 
    A2 = A.f2; 
    
    r0 = A0(varargin{:});
    r1 = A1(varargin{:});
    r2 = A2(varargin{:});

    fr = dual2(r0,r1,r2);
  else
    %fr = A(varargin{:});
    S.subs = varargin;
    fr = builtin('subsref',A,S);
  end
end

