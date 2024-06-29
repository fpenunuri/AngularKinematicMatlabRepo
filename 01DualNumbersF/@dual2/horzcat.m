%HORZCAT Horizontal concatenation for dual2 numbers.
%for Matlab

function fr = horzcat(varargin)
  vecd = builtin('horzcat',varargin{:});

  g0 = builtin('horzcat',vecd.f0);
  g1 = builtin('horzcat',vecd.f1);
  g2 = builtin('horzcat',vecd.f2);
  
  fr = dual2(g0,g1,g2);
end
