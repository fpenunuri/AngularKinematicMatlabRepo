% .*  Array multiply
%The elements of X and Y can be dual3 numbers

function local = times(varargin)
  n = nargin;    
  local = 1;
  for k = 1:n
    local = times_2args(local, varargin{k});
  end
end

function fr = times_2args(x,y)
  if (isa(x,'dual3') && isa(y,'dual3'))        
    x0 = x.f0;
    x1 = x.f1;
    x2 = x.f2;
    x3 = x.f3;
    
    y0 = y.f0;
    y1 = y.f1;
    y2 = y.f2;
    y3 = y.f3;
    
    fr0 = x0.*y0;
    fr1 = x0.*y1 + x1.*y0;
    fr2 = x0.*y2 + 2.*x1.*y1 + x2.*y0;
    fr3 = x3.*y0 + 3.*x2.*y1 + 3.*x1.*y2 + x0.*y3;
    
    fr = dual3(fr0,fr1,fr2,fr3);
    
  elseif (isa(x,'dual3') && isnumeric(y))
    x0 = x.f0;
    x1 = x.f1;
    x2 = x.f2;
    x3 = x.f3;
    
    fr0 = x0.*y;
    fr1 = x1.*y;
    fr2 = x2.*y;
    fr3 = x3.*y;
    
    fr = dual3(fr0,fr1,fr2,fr3);

  elseif (isnumeric(x) && isa(y,'dual3'))
    y0 = y.f0;
    y1 = y.f1;
    y2 = y.f2;
    y3 = y.f3;
    
    fr0 = x.*y0;
    fr1 = x.*y1;
    fr2 = x.*y2;
    fr3 = x.*y3;
    
    fr = dual3(fr0,fr1,fr2,fr3);
  else
    fr = builtin('times', x, y);
  end  
end
