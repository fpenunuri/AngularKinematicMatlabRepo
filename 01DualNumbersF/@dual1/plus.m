% X+Y; the elements of X and Y can be dual1 numbers

function local = plus (varargin)
  n = nargin;
  local = dual1(0,0);
  for k = 1:n
    local = plus_2args (local, varargin{k});
  end
end

function res = plus_2args(x,y)      
  if (isa(x,'dual1') && isa(y,'dual1'))
    
    x0 = x.f0;
    x1 = x.f1; 
    
    y0 = y.f0;
    y1 = y.f1; 
    
    r0 = builtin('plus', x0, y0);
    r1 = builtin('plus', x1, y1);
    
    res = dual1(r0,r1);
    
  elseif (isa(x,'dual1') && isnumeric(y))
    x0 = x.f0;
    x1 = x.f1; 
    
    y0 = y;
    y1 = zeros(size(y));
    
    r0 = builtin('plus', x0, y0);
    r1 = builtin('plus', x1, y1);
    
    res = dual1(r0,r1);
    
  elseif (isnumeric(x) && isa(y,'dual1')) 
    x0 = x;
    x1 = zeros(size(x));
    
    y0 = y.f0;
    y1 = y.f1; 
    
    r0 = builtin('plus', x0, y0);
    r1 = builtin('plus', x1, y1);
    
    res = dual1(r0,r1);
    
  else
    res = builtin('plus', x, y);
  end
end
