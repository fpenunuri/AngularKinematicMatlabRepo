%w = ang_vel_3pts(R1p,R0p)
%output variables
%w: angular velocity
%input variables
%R1p,R0p: matrices that store, as columns, the velocity and position, 
%respectively, of three non-collinear points.
%This function uses the ang_vel function which, in turn uses the BFF 
%method

function fr = ang_vel_3pts(R1p,R0p)
  base = @(q) base3pX(q,R1p,R0p);

  %q=t, q0p=t0=0, q1p=1=dt/dt
  fr = ang_vel(base,1,0);
end

%function that creates a basis vector as function of t
function fr = base3pX(t,R1p,R0p)
  pts = R0p + t.*R1p ;

  r1 = pts(:,1);
  r2 = pts(:,2);
  r3 = pts(:,3);

  e1 = vuni(r2-r1);
  e2 = rot_mat(pi/2,cross(e1,(r3-r1)))*e1;
  e3 = vuni(cross(e1,e2));

  fr = cat(2,e1,e2,e3);      
end

