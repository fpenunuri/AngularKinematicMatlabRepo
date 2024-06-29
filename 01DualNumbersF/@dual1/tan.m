%TAN    Tangent of argument in radians.
%The elements of X can be dual1 numbers
function fr = tan(x)
  fr = sin(x)./cos(x);
end
