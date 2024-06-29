%[w,aa] = angularKinQ12_3pts(R2p,R1p,R0p)
%output variables
%w: angular velocity
%aa: angular acceleration
%input variables
%R2p,R1p,R0p: matrices that store, as columns, the acceleration, velocity 
%and position, respectively, of three non-collinear points.
%This function uses the angularKinQ12 function which in turn uses the BFF 
%method

function [omega,alpha] = angularKinQ12_3pts(R2p,R1p,R0p)
  base = @(q) base3pX(q,R2p,R1p,R0p);

  %q2p = 0, q1p=dt/dt=1, q0p=t0=0 (q = t)
  [omega,alpha] = angularKinQ12(base,0,1,0);
end

%function that creates a basis vector as function of t
function fr = base3pX(t,R2p,R1p,R0p)
  pts = R0p + t.*R1p + 0.5*t.^2 * R2p;

  r1 = pts(:,1);
  r2 = pts(:,2);
  r3 = pts(:,3);

  e1 = vuni(r2-r1);
  e2 = rot_mat(pi/2,cross(e1,(r3-r1)))*e1;
  e3 = vuni(cross(e1,e2));

  fr = cat(2,e1,e2,e3);        
end

