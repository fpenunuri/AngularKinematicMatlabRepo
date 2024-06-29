%SIZE   Size of array.  
%for a dual2 array A, it returns size(A.f0)
function varargout = size(A,varargin)
  if(isa(A,'dual2'))
    [varargout{1:nargout}] = builtin('size',A.f0,varargin{:});
  else
    [varargout{1:nargout}] = builtin('size',A,varargin{:});     
  end
end


