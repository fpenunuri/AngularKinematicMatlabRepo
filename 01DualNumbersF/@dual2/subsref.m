% subsref function.
% A simple option that can requiere more adjustments. 
% Do not use this function inside the '@dual2' folder. For this case
% use the takepart function

function obj = subsref(obj,S)
    if(isa(obj,'dual2') && strcmp(S.type,'()'))
      obj.f0 = obj.f0(S.subs{:});
      obj.f1 = obj.f1(S.subs{:});
      obj.f2 = obj.f2(S.subs{:});
    else
      obj = builtin('subsref',obj,S);
    end
end


