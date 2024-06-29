% .*  Array multiply
%The elements of X and Y can be dual1 numbers

function local = times(varargin)
  n = nargin;    
  local = 1;
  for k = 1:n
    local = times_2args(local, varargin{k});
  end
end

function fr = times_2args(x,y)
  if (isa(x,'dual1') && isa(y,'dual1'))        
    x0 = x.f0;
    x1 = x.f1;
    
    y0 = y.f0;
    y1 = y.f1;
    
    fr0 = x0.*y0;
    fr1 = x0.*y1 + x1.*y0;
    
    fr = dual1(fr0,fr1);
    
  elseif (isa(x,'dual1') && isnumeric(y))
    x0 = x.f0;
    x1 = x.f1;
    
    fr0 = x0.*y;
    fr1 = x1.*y;
    
    fr = dual1(fr0,fr1);

  elseif (isnumeric(x) && isa(y,'dual1'))
    y0 = y.f0;
    y1 = y.f1;
    
    fr0 = x.*y0;
    fr1 = x.*y1;
    
    fr = dual1(fr0,fr1);
  else
    fr = builtin('times', x, y);
  end  
end
