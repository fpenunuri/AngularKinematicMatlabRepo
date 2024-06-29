% X+Y; the elements of X and Y can be dual2 numbers

function local = plus (varargin)
  n = nargin;
  local = dual2(0,0,0);
  for k = 1:n
    local = plus_2args (local, varargin{k});
  end
end

function res = plus_2args(x,y)      
  if (isa(x,'dual2') && isa(y,'dual2'))
    
    x0 = x.f0;
    x1 = x.f1; 
    x2 = x.f2;
    
    y0 = y.f0;
    y1 = y.f1; 
    y2 = y.f2;
    
    r0 = builtin('plus', x0, y0);
    r1 = builtin('plus', x1, y1);
    r2 = builtin('plus', x2, y2);
    
    res = dual2(r0,r1,r2);
    
  elseif (isa(x,'dual2') && isnumeric(y))
    x0 = x.f0;
    x1 = x.f1; 
    x2 = x.f2;
    
    y0 = y;
    y1 = zeros(size(y));
    y2 = y1;
    
    r0 = builtin('plus', x0, y0);
    r1 = builtin('plus', x1, y1);
    r2 = builtin('plus', x2, y2);
    
    res = dual2(r0,r1,r2);
    
  elseif (isnumeric(x) && isa(y,'dual2')) 
    x0 = x;
    x1 = zeros(size(x));
    x2 = x1;
    
    y0 = y.f0;
    y1 = y.f1; 
    y2 = y.f2;
    
    r0 = builtin('plus', x0, y0);
    r1 = builtin('plus', x1, y1);
    r2 = builtin('plus', x2, y2);
    
    res = dual2(r0,r1,r2);
    
  else
    res = builtin('plus', x, y);
  end
end
