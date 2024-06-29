%[w,aa,aj] = angularKinQ13_3pts(R3p,R2p,R1p,R0p)
%output variables
%w: angular velocity
%aa: angular acceleration
%aj: angular jerk
%input variables
%R3p,R2p,R1p,R0p: matrices that store, as columns, the jerk, acceleration, 
%velocity and position, respectively, of three non-collinear points.
%This function uses the angularKinQ13 function which in turn uses the BFF 
%method

function [omega,alpha,jerk] = ...
      angularKinQ13_3pts(R3p,R2p,R1p,R0p)

  base = @(q) base3pX(q,R3p,R2p,R1p,R0p);

  %q2p = 0, q1p=dt/dt=1, q0p=t0=0 (q = t)
  [omega,alpha,jerk] = angularKinQ13(base,0,0,1,0);
end

%function that creates a basis vector as function of t
function fr = base3pX(t,R3p,R2p,R1p,R0p)
  pts = R0p + t.*R1p + 0.5*t.^2 * R2p + t.^3 * R3p/6;

  r1 = pts(:,1);
  r2 = pts(:,2);
  r3 = pts(:,3);

  e1 = vuni(r2-r1);
  e2 = rot_mat(pi/2,cross(e1,(r3-r1)))*e1;
  e3 = vuni(cross(e1,e2));

  fr = cat(2,e1,e2,e3);       
end

