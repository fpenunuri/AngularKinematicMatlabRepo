%-  Unary minus.
%The elements of X can be dual2 numbers
function fr = uminus(x)
  if (isa(x,'dual2'))
    x0 = x.f0;
    x1 = x.f1;
    x2 = x.f2;
    
    r0 = builtin('uminus',x0);
    r1 = builtin('uminus',x1);
    r2 = builtin('uminus',x2);
    
    fr = dual2(r0,r1,r2);
  end
end 





