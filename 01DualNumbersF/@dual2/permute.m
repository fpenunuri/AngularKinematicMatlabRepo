%PERMUTE Permute array dimensions.
%   B = PERMUTE(A,ORDER) rearranges the dimensions of A so that they are in
%   the order specified by the vector ORDER.  The resulting array has the
%   same values as A but the order of the subscripts needed to access any
%   particular element is rearranged as specified by ORDER.  For an N-D
%   array A, numel(ORDER)>=ndims(A).  All the elements of ORDER must be
%   unique. A can be dual.
% 
%   PERMUTE and IPERMUTE are a generalization of transpose (.') 
%   for N-D arrays.
%
%   Example:
%      a = rand(1,2,3,4);
%      size(permute(a,[3 2 1 4])) % now it's 3-by-2-by-1-by-4.
%
%   See also SIZE  
%   IPERMUTE, CIRCSHIFT are not dualized.



function fr = permute(A,order)
A0 = A.f0;
A1 = A.f1;
A2 = A.f2;


A0d = permute(A0,order);
A1d = permute(A1,order);
A2d = permute(A2,order);

fr = dual2(A0d,A1d,A2d);
end
