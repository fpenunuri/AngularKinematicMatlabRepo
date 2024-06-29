%u = vtangent(v,w,dim) gives a unit tangent vector at 'point' v to the
%arc from v to w
%v and w are arrays of the same dimensions
%along a given dimension, elements of v and w have the same magnitude
%dim: dimension along which to operate (the same as in the sum function)
%vtangent(v,w,dim)
function fr = vtangent(varargin)
  v = varargin{1};
  w = varargin{2};

  vu = vuni(v,varargin{3:end});
  fr = vuni(w-sp(w,vu,varargin{3:end}).*vu,varargin{3:end});
end
