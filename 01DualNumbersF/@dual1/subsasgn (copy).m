% subsasgn function.
% A simple option that can requiere more adjustments.
% Do not use this function inside the '@dual1' folder. For this case
% use the setpart function 
function obj = subsasgn(obj,S,b)
    if(isa(b,'dual1') && strcmp(S.type,'()'))
        obj = dual1(obj);
        obj.f0(S.subs{:}) = b.f0;
        obj.f1(S.subs{:}) = b.f1;        
    elseif(isnumeric(b) && strcmp(S.type,'()') && isa(obj,'dual1'))
        b = dual1(b);
        obj.f0(S.subs{:}) = b.f0;
        obj.f1(S.subs{:}) = b.f1;
    else 
        obj = builtin('subsasgn',obj,S,b);
    end
end

