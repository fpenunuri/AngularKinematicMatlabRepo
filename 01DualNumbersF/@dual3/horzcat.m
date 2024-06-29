%HORZCAT Horizontal concatenation for dual3 numbers.
%for Matlab

function fr = horzcat(varargin)
  vecd = builtin('horzcat',varargin{:});

  g0 = builtin('horzcat',vecd.f0);
  g1 = builtin('horzcat',vecd.f1);
  g2 = builtin('horzcat',vecd.f2);
  g3 = builtin('horzcat',vecd.f3);
  
  fr = dual3(g0,g1,g2,g3);
end
