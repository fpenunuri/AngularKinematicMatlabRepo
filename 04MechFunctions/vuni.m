%u = vuni(v), unit vector in the direction of v
%u = vuni(v,dim), the elements of v are made into unit vectors along
%dimension dim

function fr = vuni(varargin)
  v = varargin{1};
  fr = v./sqrt(sum(v.*v,varargin{2:end}));
end
