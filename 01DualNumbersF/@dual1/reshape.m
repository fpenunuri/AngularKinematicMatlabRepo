%RESHAPE Reshape array.
%   RESHAPE(X,M,N) or RESHAPE(X,[M,N]) returns the M-by-N matrix 
%   whose elements are taken columnwise from X. An error results 
%   if X does not have M*N elements.
%
%   RESHAPE(X,M,N,P,...) or RESHAPE(X,[M,N,P,...]) returns an 
%   N-D array with the same elements as X but reshaped to have 
%   the size M-by-N-by-P-by-.... The product of the specified
%   dimensions, M*N*P*..., must be the same as NUMEL(X).
%
%   RESHAPE(X,...,[],...) calculates the length of the dimension
%   represented by [], such that the product of the dimensions 
%   equals NUMEL(X). The value of NUMEL(X) must be evenly divisible 
%   by the product of the specified dimensions. You can use only one 
%   occurrence of [].
%
%   The elements of X can be dual1 numbers

function fr = reshape(A,varargin)

  g0 = A.f0;  
  g1 = A.f1;

  fr0 = reshape(g0,varargin{:});
  fr1 = reshape(g1,varargin{:});

  fr = dual1(fr0,fr1);  
end


