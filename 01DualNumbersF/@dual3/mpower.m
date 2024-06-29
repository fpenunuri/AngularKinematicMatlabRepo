%^   Matrix power.
% A_d^x; A_d is a dual3 matrix and x an integer
% we define A_d^0 as a diagonal matrix with (A_d)ii=dual3(1,0,0,0) 
function fr = mpower(A,x)
  if(isa(A,'dual3') && ~isa(x,'dual3') && length(A)*length(x) ~=1)
    if(isa(A,'dual3') && length(x)==1  && x == 0)
      A0 = A.f0;
      uno = eye(size(A0));
      cero = zeros(size(A0));
      fr = dual3(uno,cero,cero,cero);
      
    elseif(isa(A,'dual3') && length(x) == 1 && mod(x,1) == 0 && x > 0)
      fr = A;
      for k = 1:(x-1)
        fr = fr*A;
      end
      
    elseif(isa(A,'dual3') && length(x) == 1 && mod(x,1) == 0 && x < 0 )
      fr = A;
      for k = 1:(-x - 1)
        fr = fr*A;
      end
      fr = inv(fr);

    elseif(isa(A,'dual3') && (mod(x,1) ~= 0 || length(x) ~= 1) )
      error('dual3 mpower is only supported for integer exponents.') 
    end
    
  elseif(isa(A,'dual3') || isa(x,'dual3') && length(x)*length(A)==1)
    fr = A.^x; 
  end
end
