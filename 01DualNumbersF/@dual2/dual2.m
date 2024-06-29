%Dual2 number class

%F. Penunuri
%UADY
classdef dual2
%Class for dual numbers 
  properties
    f0;%function
    f1;%first derivative
    f2;%second derivative
  end
  
  methods
    function obj = dual2(g0,g1,g2)
      if (nargin == 1)
        if (isa(g0,'dual2'))
          obj=g0;                   
        elseif(isnumeric(g0))
          obj.f0 = g0;
          obj.f1 = zeros(size(g0));
          obj.f2 = zeros(size(g0));
        else
          error(['when only one argument is supplied to the dual2'...
                 ' constructor, this must be numeric or dual2'])
        end
      elseif(nargin == 2)
        error('wrong number of arguments in dual2 constructor')          
      elseif(nargin == 3)               
        if(isnumeric(g0) && isequal(size(g0),size(g1),size(g2)))
          obj.f0 = g0;
          obj.f1 = g1;
          obj.f2 = g2;                    
        else
          error(['arguments in dual2 constructor must be numeric'...
                 ' having the same dimensions'])
        end
      elseif(nargin > 3)
        error('too many input arguments in dual2 constructor')
      else
        error('wrong number of arguments in dual2 constructor')
      end
    end
  end
  
  methods (Hidden)
    function z = zerosLike(obj,varargin)
      if(isa(obj,'dual2'))
        z = dual2(zeros(varargin{:}));
      end
    end
  end
end
