% X+Y; the elements of X and Y can be dual3 numbers

function local = plus (varargin)
  n = nargin;
  local = dual3(0,0,0,0);
  for k = 1:n
    local = plus_2args (local, varargin{k});
  end
end

function res = plus_2args(x,y)      
  if (isa(x,'dual3') && isa(y,'dual3'))
    
    x0 = x.f0;
    x1 = x.f1; 
    x2 = x.f2;
    x3 = x.f3;
    
    y0 = y.f0;
    y1 = y.f1; 
    y2 = y.f2;
    y3 = y.f3;
    
    r0 = builtin('plus', x0, y0);
    r1 = builtin('plus', x1, y1);
    r2 = builtin('plus', x2, y2);
    r3 = builtin('plus', x3, y3);
    
    res = dual3(r0,r1,r2,r3);
    
  elseif (isa(x,'dual3') && isnumeric(y))
    x0 = x.f0;
    x1 = x.f1; 
    x2 = x.f2;
    x3 = x.f3;
    
    y0 = y;
    y1 = zeros(size(y));
    y2 = y1;
    y3 = y1;
    
    r0 = builtin('plus', x0, y0);
    r1 = builtin('plus', x1, y1);
    r2 = builtin('plus', x2, y2);
    r3 = builtin('plus', x3, y3);
    
    res = dual3(r0,r1,r2,r3);
    
  elseif (isnumeric(x) && isa(y,'dual3')) 
    x0 = x;
    x1 = zeros(size(x));
    x2 = x1;
    x3 = x1;
    
    y0 = y.f0;
    y1 = y.f1; 
    y2 = y.f2;
    y3 = y.f3;
    
    r0 = builtin('plus', x0, y0);
    r1 = builtin('plus', x1, y1);
    r2 = builtin('plus', x2, y2);
    r3 = builtin('plus', x3, y3);
    
    res = dual3(r0,r1,r2,r3);
    
  else
    res = builtin('plus', x, y);
  end
end
