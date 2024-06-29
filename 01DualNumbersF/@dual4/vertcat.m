%VERTCAT vertical concatenation for dual4 numbers
% For Matlab

function fr = vertcat(varargin)
  vecd = builtin('vertcat',varargin{:});
  
  g0 = builtin('vertcat',vecd.f0);
  g1 = builtin('vertcat',vecd.f1);
  g2 = builtin('vertcat',vecd.f2);
  g3 = builtin('vertcat',vecd.f3);
  g4 = builtin('vertcat',vecd.f4);
  
  fr = dual4(g0,g1,g2,g3,g4);
end

