%^   Matrix power.
% A_d^x; A_d is a dual1 matrix and x an integer
% we define A_d^0 as a diagonal matrix with (A_d)ii=dual1(1,0) 
function fr = mpower(A,x)
  if(isa(A,'dual1') && ~isa(x,'dual1') && length(A)*length(x) ~=1)
    if(isa(A,'dual1') && length(x)==1  && x == 0)
      A0 = A.f0;
      uno = eye(size(A0));
      cero = zeros(size(A0));
      fr = dual1(uno,cero);
      
    elseif(isa(A,'dual1') && length(x) == 1 && mod(x,1) == 0 && x > 0)
      fr = A;
      for k = 1:(x-1)
        fr = fr*A;
      end
      
    elseif(isa(A,'dual1') && length(x) == 1 && mod(x,1) == 0 && x < 0 )
      fr = A;
      for k = 1:(-x - 1)
        fr = fr*A;
      end
      fr = inv(fr);

    elseif(isa(A,'dual1') && (mod(x,1) ~= 0 || length(x) ~= 1) )
      error('dual1 mpower is only supported for integer exponents.') 
    end
    
  elseif(isa(A,'dual1') || isa(x,'dual1') && length(x)*length(A)==1)
      fr = A.^x;
  else
      fr = builtin('mpower', A, x);
  end

end
