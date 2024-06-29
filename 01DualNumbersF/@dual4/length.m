%LENGTH   Length of vector.
%   LENGTH(X) returns the length of vector X.  It is equivalent
%   to MAX(SIZE(X)) for non-empty arrays and 0 for empty ones.
%  
% For a dual4 array A, it returns length(A.f0)
function fr = length(A)
  if(isa(A,'dual4'))
    fr = builtin('length',A.f0);
  else
    fr = builtin('length',A);
  end
end


