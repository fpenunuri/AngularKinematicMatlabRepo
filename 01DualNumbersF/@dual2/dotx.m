%dot(X,Y,dim) product
%The elements of X and Y can be dual2 numbers

%F. Penunuri
%UADY
function fr = dotx(varargin)
  [args{1:nargin}] = varargin{:};
  arr1 = args{1};
  arr2 = args{2};

  A0 = arr1.f0; A1 = arr1.f1; A2 = arr1.f2;
  B0 = arr2.f0; B1 = arr2.f1; B2 = arr2.f2;

%   A0 = zeros(2,3); A1=A0; A2=A0;
%   B0 = zeros(2,3); B1=B0; B2=B0;
  
  
  args0 = args;
  args0{1} = A0; args0{2} = B0;

  args10 = args;
  args10{1} = A1; args10{2} = B0;
  
  args01 = args;
  args01{1} = A0; args01{2} = B1;
  
  args02 = args;
  args02{1} = A0; args02{2} = B2;
  
  args11 = args;
  args11{1} = A1; args11{2} = B1;

  args20 = args;
  args20{1} = A2; args20{2} = B0;  
  
  fr0 = dot(args0{:});
  fr1 = dot(args10{:}) + dot(args01{:});
  fr2 = dot(args02{:}) + 2.*dot(args11{:}) + dot(args20{:});

%   fr0 = builtin('dot',args0{:});
%   fr1 = builtin('dot',args10{:}) + builtin('dot',args01{:});
%   fr2 = builtin('dot',args02{:}) + 2.*builtin('dot',args11{:}) + ...
%       builtin('dot',args20{:});
   
  fr = dual2(fr0,fr1,fr2);
end
