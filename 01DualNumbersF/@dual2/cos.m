%COS    Cosine of argument in radians.
% The elements of X can be dual2 numbers
function fr = cos(x)
% fr = sin(x + pi/2);
        
    g0 = x.f0;
    g1 = x.f1;
    g2 = x.f2;
    
    cos0 = cos(g0);
    cos1 = -sin(g0).*g1;
    cos2 = -cos(g0).*(g1.^2) - sin(g0).*g2;
    
    fr = dual2(cos0,cos1,cos2);	
end
