%ABSX: Absolute value without conjugation. Intended to be used with the
%complex steep approximation method
%   ABSX(X)
% The elements of X can be dual1 numbers. 
function fr = absX(g)
  fr = sqrt(g.*g);
end
