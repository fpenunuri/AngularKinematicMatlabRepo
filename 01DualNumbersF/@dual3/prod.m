%PROD Product of elements.
%the elements of X can be dual3 numbers
%for the dual case only the prod(X,dim) is implemented
%for the product of all the elements of X use prod(X(:),1)
%or use the prod_elements function
%This function could be very time consuming for matrices of rank (number 
%of indices) greater than 2

function local = prod(A, dim)
  A0 = A.f0;
  sA  = size(A0);
  if (isvector(A0)) %1*N or N*1 array
    if(sA(dim) == max(sA))
      local = prod_elements(A);
    elseif(mod(dim,1)==0 && dim > 0)
      local = A;
    else
      error('invalid dimension')
    end
  elseif(isvector(squeeze(A0)))
    if(sA(dim) == max(sA))
      local = prod_elements(A);
    elseif(mod(dim,1)==0 && dim > 0)
      local = A;
    else  %esta parte nunca se alcanza, solo esta por completez
      error('invalid dimension')
    end
  else
    Ra = length(sA);
    indx = cell(1,Ra);
    [indx{:}] = deal(':');        
    local = 1;
    for k = 1:sA(dim)
      indx{dim} = k;
      local = local.*takepart(A,indx{:});
    end
  end
end

