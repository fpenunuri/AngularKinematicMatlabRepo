%-  Unary minus.
%The elements of X can be dual1 numbers
function fr = uminus(x)
  if (isa(x,'dual1'))
    x0 = x.f0;
    x1 = x.f1;
    
    r0 = builtin('uminus',x0);
    r1 = builtin('uminus',x1);
    
    fr = dual1(r0,r1);
  end
end 





