%Atan2 for complex arguments

function fr = atan2_z(y,x)
    fr = zeros(size(y));
    
    log1 = x >  0;
    log2 = x <  0 & y >= 0;
    log3 = x <  0 & y < 0;
    log4 = x == 0 & y > 0;
    log5 = x == 0 & y < 0;
    log6 = x == 0 & y == 0;
    
    fr(log6) = NaN;
    fr(log5) = -pi/2;
    fr(log4) = pi/2;
 
    A1 = atan(y./x);
    A2 = atan(y./x) + pi;
    A3 = atan(y./x) - pi;
    fr(log1) = A1(log1);
    fr(log2) = A2(log2);
    fr(log3) = A3(log3);
end

% $$$ although correct, this version present precision problems 
% $$$   
% $$$   function  f_res = atan2_z(zy,zx) 
% $$$     ii=1i;
% $$$     
% $$$     x1 = real(zx);
% $$$     x2 = imag(zx);
% $$$ 
% $$$     y1 = real(zy);
% $$$     y2 = imag(zy);
% $$$ 
% $$$     num = x1 + ii*x2 + ii*(y1 + ii*y2);
% $$$     den = sqrt((x1 + ii*x2).^2 + (y1 + ii*y2).^2);
% $$$     divnd = num./den;
% $$$     t1 = atan2(imag(divnd),real(divnd));
% $$$     t2 = ii*log(sqrt((x2 + y1).^2 + (x1 - y2).^2)./((2*x1.*x2 +  ...
% $$$          2*y1.*y2).^2 + (x1.^2 - x2.^2 + y1.^2 - y2.^2).^2).^(1/4));
% $$$ 
% $$$     f_res = t1 - t2;
% $$$   end
