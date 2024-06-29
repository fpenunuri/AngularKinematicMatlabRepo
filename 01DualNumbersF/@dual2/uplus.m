%+  Unary plus.
%The elements of X can be dual2 numbers
function fr = uplus(x)
  if (isa(x,'dual2'))
    fr = x;
  end
end 