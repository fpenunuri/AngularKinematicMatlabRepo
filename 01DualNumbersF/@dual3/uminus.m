%-  Unary minus.
%The elements of X can be dual3 numbers
function fr = uminus(x)
  if (isa(x,'dual3'))
    x0 = x.f0;
    x1 = x.f1;
    x2 = x.f2;
    x3 = x.f3;
    
    r0 = builtin('uminus',x0);
    r1 = builtin('uminus',x1);
    r2 = builtin('uminus',x2);
    r3 = builtin('uminus',x3);
    
    fr = dual3(r0,r1,r2,r3);
  end
end 





