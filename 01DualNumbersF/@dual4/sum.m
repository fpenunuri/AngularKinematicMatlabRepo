%SUM Sum of elements.
%   S = SUM(X) is the sum of the elements of the vector X. If X is a matrix,
%   S is a row vector with the sum over each column. For N-D arrays, 
%   SUM(X) operates along the first non-singleton dimension.
%
%   S = SUM(X,DIM) sums along the dimension DIM. 
%The elements of X can be dual4 numbers

%F. Penunuri
%UADY
function fr = sum(varargin)
  [args{1:nargin}] = varargin{:};
  arr = args{1};

  g0 = arr.f0;  
  g1 = arr.f1;
  g2 = arr.f2;
  g3 = arr.f3;
  g4 = arr.f4;

  args0 = args;
  args0{1} = g0;

  args1 = args;
  args1{1} = g1;

  args2 = args;
  args2{1} = g2;

  args3 = args;
  args3{1} = g3;

  args4 = args;
  args4{1} = g4;
  
  fr0 = sum(args0{:});
  fr1 = sum(args1{:});
  fr2 = sum(args2{:});
  fr3 = sum(args3{:});
  fr4 = sum(args4{:});

  fr = dual4(fr0,fr1,fr2,fr3,fr4);
end






