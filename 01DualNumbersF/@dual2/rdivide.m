%./  Right array divide.
% The elements of X and Y can be dual2 numbers
function fr = rdivide(x,y)
  if (isa(x,'dual2') && isa(y,'dual2'))
    x0 = x.f0;
    x1 = x.f1;
    x2 = x.f2;
    
    y0 = y.f0;
    y1 = y.f1;
    y2 = y.f2;
    
    
    fr0 = x0./y0;
    fr1 = (x1.*y0 - x0.*y1)./y0.^2;
    fr2 = (x2.*y0.^2 - 2*x1.*y0.*y1 + 2*x0.*y1.^2 - x0.*y0.* ...
           y2)./(y0.^3);
    
    fr = dual2(fr0,fr1,fr2);
    
  elseif (isa(x,'dual2') && isnumeric(y))
    x0 = x.f0;
    x1 = x.f1;
    x2 = x.f2;
    
    fr0 = x0./y; 
    fr1 = x1./y;
    fr2 = x2./y;
    
    fr = dual2(fr0,fr1,fr2);
    
  elseif (isnumeric(x) && isa(y,'dual2'))
    y0 = y.f0;
    y1 = y.f1;
    y2 = y.f2;
    
    fr0 = x./y0;
    fr1 = (-y1.*x)./(y0.^2);
    fr2 = x.*(2.*y1.^2 - y0.*y2)./(y0.^3);
    
    fr = dual2(fr0,fr1,fr2);
  end
end  
