%+  Unary plus.
%The elements of X can be dual numbers
function fr = uplus(x)
  if (isa(x,'dual1'))
    fr = x;
  end
end 