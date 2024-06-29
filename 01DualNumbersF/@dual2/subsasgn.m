% subsasgn function.
% A simple option that can requiere more adjustments.
% Do not use this function inside the '@dual2' folder. For this case
% use the setpart function 
function obj = subsasgn(obj,S,b)
    if(isa(b,'dual2') && strcmp(S.type,'()'))
        obj = dual2(obj);
        obj.f0(S.subs{:}) = b.f0;
        obj.f1(S.subs{:}) = b.f1;
        obj.f2(S.subs{:}) = b.f2;        
    elseif(isnumeric(b) && strcmp(S.type,'()') && isa(obj,'dual2'))
        b = dual2(b);
        obj.f0(S.subs{:}) = b.f0;
        obj.f1(S.subs{:}) = b.f1;
        obj.f2(S.subs{:}) = b.f2;
    else 
        obj = builtin('subsasgn',obj,S,b);
    end
end

