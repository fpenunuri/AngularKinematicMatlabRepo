%Dual3 number class

%F. Penunuri
%UADY
classdef dual3
%Class for dual numbers 
  properties
    f0;
    f1;
    f2;
    f3;
  end
  
  methods
    function obj = dual3(g0,g1,g2,g3)
      if (nargin == 1)
        if (isa(g0,'dual3'))
          obj=g0;                   
        elseif(isnumeric(g0))
          obj.f0 = g0;
          obj.f1 = zeros(size(g0));
          obj.f2 = zeros(size(g0));
          obj.f3 = zeros(size(g0));
        else
          error(['when only one argument is supplied to the dual3'...
                 ' constructor, this must be numeric or dual3'])
        end
      elseif(nargin == 2)
        error('wrong number of arguments in dual3 constructor')

      elseif(nargin == 3)
        error('wrong number of arguments in dual3 constructor')
        
      elseif(nargin == 4)               
        if(isnumeric(g0) && isequal(size(g0),size(g1),size(g2),...
                                    size(g3)))
          obj.f0 = g0;
          obj.f1 = g1;
          obj.f2 = g2;
          obj.f3 = g3;
        else
          error(['arguments in dual3 constructor must be numeric'...
                 ' having the same dimensions'])
        end
      elseif(nargin > 4)
        error('too many input arguments in dual3 constructor')
      else
        error('wrong number of arguments in dual3 constructor')
      end
    end
  end
  
  methods (Hidden)
    function z = zerosLike(obj,varargin)
      if(isa(obj,'dual3'))
        z = dual3(zeros(varargin{:}));
      end
    end
  end    
end
