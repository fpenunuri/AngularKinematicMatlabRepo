%Dual4 number class

%F. Penunuri
%UADY
classdef dual4
%Class for dual4 numbers 
  properties
    f0; f1; f2; f3; f4;
  end
  
  methods
    function obj = dual4(g0,g1,g2,g3,g4)
      if (nargin == 1)
        if (isa(g0,'dual4'))
          obj=g0;                   
        elseif(isnumeric(g0))
          obj.f0 = g0;
          obj.f1 = zeros(size(g0));
          obj.f2 = zeros(size(g0));
          obj.f3 = zeros(size(g0));
          obj.f4 = zeros(size(g0));
        else
          error('operation not supported in dual4 constructor')
        end
      elseif(nargin == 2)
        error('wrong number of arguments in dual4 constructor')

      elseif(nargin == 3)
        error('wrong number of arguments in dual4 constructor')

      elseif(nargin == 4)
        error('wrong number of arguments in dual4 constructor')
        
      elseif(nargin == 5)               
        if(isnumeric(g0) && isequal(size(g0),size(g1),size(g2),...
                                    size(g3),size(g4)))
          obj.f0 = g0;
          obj.f1 = g1;
          obj.f2 = g2;
          obj.f3 = g3;
          obj.f4 = g4;
        else
          error(['arguments in dual4 constructor must be numeric'...
                 ' having the same dimensions'])
        end
      elseif(nargin > 5)
        error('too many input arguments in dual4 constructor')
      else
        error('wrong number of arguments in dual4 constructor')
      end
    end
  end
  methods (Hidden)
    function z = zerosLike(obj,varargin)
      if(isa(obj,'dual4'))
        z = dual4(zeros(varargin{:}));
      end
    end
  end    
end
