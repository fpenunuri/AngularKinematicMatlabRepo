%+  Unary plus.
%The elements of X can be dual4 numbers

function fr = uplus(x)
  if (isa(x,'dual4'))
    fr = x;
  end
end
