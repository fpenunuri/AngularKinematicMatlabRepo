%+  Unary plus.
%The elements of X can be dual3 numbers
function fr = uplus(x)
  if (isa(x,'dual3'))
    fr = x;
  end
end
