%fd = f2dualf(f)
%
%If any of the arguments of f are dual, this function will convert f into
%a function fd with dual arguments.

function frv = f2dualf(f)
  frv = @(varargin) faux(f,varargin{:});  
end

function fr = faux(ff,varargin)
  vars = args2duals(varargin{:});

  fr = ff(vars{:});
end
