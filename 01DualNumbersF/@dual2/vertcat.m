%VERTCAT vertical concatenation for dual2 numbers
% For Matlab

function fr = vertcat(varargin)
  vecd = builtin('vertcat',varargin{:});
  
  g0 = builtin('vertcat',vecd.f0);
  g1 = builtin('vertcat',vecd.f1);
  g2 = builtin('vertcat',vecd.f2);
  
  fr = dual2(g0,g1,g2);
end

