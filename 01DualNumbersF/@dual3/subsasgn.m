% subsasgn function.
% A simple option that can requiere more adjustments.
% Do not use this function inside the '@dual3' folder. For this case
% use the setpart function 
function obj = subsasgn(obj,S,b)
    if(isa(b,'dual3') && strcmp(S.type,'()'))
        obj = dual3(obj);
        obj.f0(S.subs{:}) = b.f0;
        obj.f1(S.subs{:}) = b.f1;
        obj.f2(S.subs{:}) = b.f2;        
        obj.f3(S.subs{:}) = b.f3;        
    elseif(isnumeric(b) && strcmp(S.type,'()') && isa(obj,'dual3'))
        b = dual3(b);
        obj.f0(S.subs{:}) = b.f0;
        obj.f1(S.subs{:}) = b.f1;
        obj.f2(S.subs{:}) = b.f2;
        obj.f3(S.subs{:}) = b.f3;
    else 
        obj = builtin('subsasgn',obj,S,b);
    end
end

