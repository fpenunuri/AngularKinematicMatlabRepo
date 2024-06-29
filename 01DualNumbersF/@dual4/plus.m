% X+Y; the elements of X and Y can be dual4 numbers

function local = plus (varargin)
  n = nargin;
  local = dual4(0,0,0,0,0);
  for k = 1:n
    local = plus_2args (local, varargin{k});
  end
end

function res = plus_2args(x,y)      
  if (isa(x,'dual4') && isa(y,'dual4'))
    
    x0 = x.f0;
    x1 = x.f1; 
    x2 = x.f2;
    x3 = x.f3;
    x4 = x.f4;
    
    y0 = y.f0;
    y1 = y.f1; 
    y2 = y.f2;
    y3 = y.f3;
    y4 = y.f4;
    
    r0 = builtin('plus', x0, y0);
    r1 = builtin('plus', x1, y1);
    r2 = builtin('plus', x2, y2);
    r3 = builtin('plus', x3, y3);
    r4 = builtin('plus', x4, y4);
    
    res = dual4(r0,r1,r2,r3,r4);
    
  elseif (isa(x,'dual4') && isnumeric(y))
    x0 = x.f0;
    x1 = x.f1; 
    x2 = x.f2;
    x3 = x.f3;
    x4 = x.f4;
    
    y0 = y;
    y1 = zeros(size(y));
    y2 = y1;
    y3 = y1;
    y4 = y1;
    
    r0 = builtin('plus', x0, y0);
    r1 = builtin('plus', x1, y1);
    r2 = builtin('plus', x2, y2);
    r3 = builtin('plus', x3, y3);
    r4 = builtin('plus', x4, y4);
    
    res = dual4(r0,r1,r2,r3,r4);
    
  elseif (isnumeric(x) && isa(y,'dual4')) 
    x0 = x;
    x1 = zeros(size(x));
    x2 = x1;
    x3 = x1;
    x4 = x1;
    
    y0 = y.f0;
    y1 = y.f1; 
    y2 = y.f2;
    y3 = y.f3;
    y4 = y.f4;
    
    r0 = builtin('plus', x0, y0);
    r1 = builtin('plus', x1, y1);
    r2 = builtin('plus', x2, y2);
    r3 = builtin('plus', x3, y3);
    r4 = builtin('plus', x4, y4);
    
    res = dual4(r0,r1,r2,r3,r4);
    
  else
    res = builtin('plus', x, y);
  end
end
