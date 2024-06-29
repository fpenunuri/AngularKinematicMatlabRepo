%HORZCAT Horizontal concatenation for dual4 numbers.
%for Matlab

function fr = horzcat(varargin)
  vecd = builtin('horzcat',varargin{:});

  g0 = builtin('horzcat',vecd.f0);
  g1 = builtin('horzcat',vecd.f1);
  g2 = builtin('horzcat',vecd.f2);
  g3 = builtin('horzcat',vecd.f3);
  g4 = builtin('horzcat',vecd.f4);
  
  fr = dual4(g0,g1,g2,g3,g4);
end
