% subsasgn function.
% A simple option that can requiere more adjustments.
% Do not use this function inside the '@dual4' folder. For this case
% use the setpart function

function obj = subsasgn(obj,S,b)
    if(isa(b,'dual4') && strcmp(S.type,'()'))
        obj = dual4(obj);
        obj.f0(S.subs{:}) = b.f0;
        obj.f1(S.subs{:}) = b.f1;
        obj.f2(S.subs{:}) = b.f2;        
        obj.f3(S.subs{:}) = b.f3;
        obj.f4(S.subs{:}) = b.f4;
    elseif(isnumeric(b) && strcmp(S.type,'()') && isa(obj,'dual4'))
        b = dual4(b);
        obj.f0(S.subs{:}) = b.f0;
        obj.f1(S.subs{:}) = b.f1;
        obj.f2(S.subs{:}) = b.f2;
        obj.f3(S.subs{:}) = b.f3;
        obj.f4(S.subs{:}) = b.f4;
    else 
        obj = builtin('subsasgn',obj,S,b);
    end
end


