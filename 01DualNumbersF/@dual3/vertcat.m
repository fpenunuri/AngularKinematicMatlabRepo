%VERTCAT vertical concatenation for dual3 numbers
% For Matlab

function fr = vertcat(varargin)
  vecd = builtin('vertcat',varargin{:});
  
  g0 = builtin('vertcat',vecd.f0);
  g1 = builtin('vertcat',vecd.f1);
  g2 = builtin('vertcat',vecd.f2);
  g3 = builtin('vertcat',vecd.f3);
  
  fr = dual3(g0,g1,g2,g3);
end

