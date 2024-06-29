%use this function to select elements, from a dual1 number 
%with arrays components
%Example: takepart(A,1,2) = A(1,2)
function fr = takepart(A,varargin)
  if(isa(A,'dual1'))        
    A0 = A.f0;
    A1 = A.f1; 
    
    r0 = A0(varargin{:});
    r1 = A1(varargin{:});

    fr = dual1(r0,r1);
  else
    S.subs = varargin;
    fr = builtin('subsref',A,S);
  end
end

