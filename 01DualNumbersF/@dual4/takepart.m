%use this function to select elements, from a dual4 number 
%with arrays components
%Example: takepart(A,1,2) = A(1,2)

function fr = takepart(A,varargin)
  if(isa(A,'dual4'))        
    A0 = A.f0;
    A1 = A.f1; 
    A2 = A.f2; 
    A3 = A.f3; 
    A4 = A.f4; 
    
    r0 = A0(varargin{:});
    r1 = A1(varargin{:});
    r2 = A2(varargin{:});
    r3 = A3(varargin{:});
    r4 = A4(varargin{:});

    fr = dual4(r0,r1,r2,r3,r4);
  else
    S.subs = varargin;
    fr = builtin('subsref',A,S);
  end
end

