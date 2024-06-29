%VERTCAT vertical concatenation for dual1 numbers
% For Matlab

function fr = vertcat(varargin)
  vecd = builtin('vertcat',varargin{:});
  
  g0 = builtin('vertcat',vecd.f0);
  g1 = builtin('vertcat',vecd.f1);
  
  fr = dual1(g0,g1);
end

