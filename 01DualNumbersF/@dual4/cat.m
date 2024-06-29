%CAT Concatenate arrays.
%   CAT(DIM,A,B) concatenates the arrays A and B along
%   the dimension DIM. 
% The arrays can be dual4 numbers

%for Matlab

function fr = cat(n,varargin)
  vecd = builtin('horzcat',varargin{:});

  g0 = {vecd.f0};
  g1 = {vecd.f1};
  g2 = {vecd.f2};
  g3 = {vecd.f3};
  g4 = {vecd.f4};
  
  g0 =  builtin('cat',n,g0{:});
  g1 =  builtin('cat',n,g1{:});
  g2 =  builtin('cat',n,g2{:});
  g3 =  builtin('cat',n,g3{:});
  g4 =  builtin('cat',n,g4{:});
  
  fr = dual4(g0,g1,g2,g3,g4);
end
