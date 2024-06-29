%COS    Cosine of argument in radians.
% The elements of X can be dual1 numbers
function fr = cos(x)
% fr = sin(x + pi/2);
    
    g0 = x.f0;
    g1 = x.f1;
    
    cos0 = cos(g0);
    cos1 = -sin(g0).*g1;
    
    fr = dual1(cos0,cos1);	
end
