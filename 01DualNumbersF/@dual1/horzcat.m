%HORZCAT Horizontal concatenation for dual1 numbers.
%for Matlab

function fr = horzcat(varargin)
  vecd = builtin('horzcat',varargin{:});

  g0 = builtin('horzcat',vecd.f0);
  g1 = builtin('horzcat',vecd.f1);
  
  fr = dual1(g0,g1);
end
