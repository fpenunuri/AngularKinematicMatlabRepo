%-  Unary minus.
%The elements of X can be dual4 numbers

function fr = uminus(x)
  if (isa(x,'dual4'))
    x0 = x.f0;
    x1 = x.f1;
    x2 = x.f2;
    x3 = x.f3;
    x4 = x.f4;
    
    r0 = builtin('uminus',x0);
    r1 = builtin('uminus',x1);
    r2 = builtin('uminus',x2);
    r3 = builtin('uminus',x3);
    r4 = builtin('uminus',x4);
    
    fr = dual4(r0,r1,r2,r3,r4);
  end
end 





