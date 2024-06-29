addpath('./01DualNumbersF;./02AuxiliarDualFun;./03KinematicQuantities')
addpath('./04MechFunctions')

%Calculation of velocity and acceleration from three non-collinear points
%for two classical examples

%EXAMPLE A
p10p = [100;100;0]; 
p20p = [300;300;0]; 
p30p = [220;180;0]; 

p11p = [600;-400;100]; 
p21p = [200;0;0]; 
p31p = [440;-160;40]; 

p12p = [850;1200;-240]; 
p22p = [200;200;0]; 
p32p = [420;760;-140];

R0p = cat(2,p10p,p20p,p30p);
R1p = cat(2,p11p,p21p,p31p);
R2p = cat(2,p12p,p22p,p32p);

w = ang_vel_3pts_vec(R1p,R0p);
a = ang_accel_3pts_vec(R2p,R1p,R0p);  

fprintf("Angular kinematic quantities (AKQ's) as columns\n\n")

disp("AKQ's for Example 1")
fprintf('%s\t%s\n', 'velocity', 'acceleration')
disp([w,a])

%--------------------------------------------------------------------------
%EXAMPLE B
p10p = [1/2; -sqrt(3)/6; 0]; 
p20p = [0; sqrt(3)/3; 0]; 
p30p = [-1/2; -sqrt(3)/6; 0]; 

p11p = (4 - sqrt(2))/4 * [0;0;1]; 
p21p = (4 - sqrt(3))/4 * [0;0;1]; 
p31p = (4 + sqrt(2))/4 * [0;0;1]; 

p12p = [-6 + 4*sqrt(3); 12 - 3*sqrt(2); 0]/24; 
p22p = -[8*sqrt(3) + 3*sqrt(6); 3*sqrt(3); 0]/24; 
p32p = [6 + 4*sqrt(3); -12 + 3*sqrt(2); 0]/24; 

%a = p31p - p11p; c = p20p - p10p;
%str = sprintf('%s%0.5f','a.c = ', a.' * c);
%disp(str)

R0p = cat(2,p10p,p20p,p30p);
R1p = cat(2,p11p,p21p,p31p);
R2p = cat(2,p12p,p22p,p32p);

disp("AKQ's for Example 2")
fprintf('%s\t%s\n', 'velocity', 'acceleration')
[w,a] = angularKinQ12_3pts(R2p,R1p,R0p);
disp([w,a])

