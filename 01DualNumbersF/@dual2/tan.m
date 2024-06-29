%TAN    Tangent of argument in radians.
%The elements of X can be dual2 numbers
function fr = tan(x)
  fr = sin(x)./cos(x);
end
