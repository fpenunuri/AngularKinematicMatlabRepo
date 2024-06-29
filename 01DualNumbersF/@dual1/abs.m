%ABS    Absolute value.
%   ABS(X)
% The elements of X can be dual1 numbers. 
function fr = abs(g)
    fr = sqrt(g.*conj(g));
end
