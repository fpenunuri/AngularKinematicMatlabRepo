% .*  Array multiply
%The elements of X and Y can be dual4 numbers


function local = times(varargin)
  n = nargin;    
  local = 1;
  for k = 1:n
    local = times_2args(local, varargin{k});
  end
end

function fr = times_2args(x,y)
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
    
    fr0 = x0.*y0;
    fr1 = x0.*y1 + x1.*y0;
    fr2 = x0.*y2 + 2.*x1.*y1 + x2.*y0;
    fr3 = x3.*y0 + 3.*x2.*y1 + 3.*x1.*y2 + x0.*y3;
    fr4 = x4.*y0 + 4.*x3.*y1 + 6.*x2.*y2 + 4.*x1.*y3 + x0.*y4;
    
    fr = dual4(fr0,fr1,fr2,fr3,fr4);
    
  elseif (isa(x,'dual4') && isnumeric(y))
    x0 = x.f0;
    x1 = x.f1;
    x2 = x.f2;
    x3 = x.f3;
    x4 = x.f4;
    
    fr0 = x0.*y;
    fr1 = x1.*y;
    fr2 = x2.*y;
    fr3 = x3.*y;
    fr4 = x4.*y;
    
    fr = dual4(fr0,fr1,fr2,fr3,fr4);

  elseif (isnumeric(x) && isa(y,'dual4'))
    y0 = y.f0;
    y1 = y.f1;
    y2 = y.f2;
    y3 = y.f3;
    y4 = y.f4;
    
    fr0 = x.*y0;
    fr1 = x.*y1;
    fr2 = x.*y2;
    fr3 = x.*y3;
    fr4 = x.*y4;
    
    fr = dual4(fr0,fr1,fr2,fr3,fr4);
  else
    fr = builtin('times', x, y);
  end  
end
