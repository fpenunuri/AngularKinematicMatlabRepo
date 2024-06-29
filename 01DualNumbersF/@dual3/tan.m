%TAN    Tangent of argument in radians.
%The elements of X can be dual3 numbers
function fr = tan(x)
  fr = sin(x)./cos(x);
end
