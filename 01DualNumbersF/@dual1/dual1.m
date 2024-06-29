%Dual1 number class

%F. Penunuri
%UADY
%Mexico, July 2022
classdef dual1
%Class for dual numbers 
  properties
    f0;%function
    f1;%first derivative
  end
  
  methods
    function obj = dual1(g0,g1)
      if (nargin == 1)
        if (isa(g0,'dual1'))
          obj=g0;                   
        elseif(isnumeric(g0))
          obj.f0 = g0;
          obj.f1 = zeros(size(g0));
        else
          error('Operation not supported in dual1 constructor')
        end
      elseif(nargin == 2)               
        if(isnumeric(g0) && isequal(size(g0),size(g1)))
          obj.f0 = g0;
          obj.f1 = g1;
        else
          error(['arguments in dual1 constructor must be numeric'...
                 ' having the same dimensions'])
        end
      elseif(nargin > 2)
        error('too many input arguments in dual1 constructor')
      else
        error('wrong number of arguments in dual1 constructor')
      end
    end       
  end
  
  methods (Hidden)
    function z = zerosLike(obj,varargin)
      if(isa(obj,'dual1'))
        z = dual1(zeros(varargin{:}));
      end
    end
  end
end
