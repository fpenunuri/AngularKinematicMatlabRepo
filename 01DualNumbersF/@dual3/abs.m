%ABS    Absolute value.
%   ABS(X)
% The elements of X can be dual3 numbers. 
function fr = abs(g)
    fr = sqrt(g.*conj(g));
end
