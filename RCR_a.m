addpath('./01DualNumbersF;./02AuxiliarDualFun;./03KinematicQuantities')
addpath('./04MechFunctions')

%--------------------------------------------------------------------------
%The RCR robot manipulator
%Angular Kinematic quantities for a given set of generalized coordinates

th0p = pi/2; phi0p = 0; s0p = 2; bt0p = 0;
th1p = 1;    phi1p = 5; s1p = 1; bt1p = 1;
th2p = 1;    phi2p = 0; s2p = 2; bt2p = 1;
th3p = 1;    phi3p = 2; s3p = 3; bt3p = 4;
th4p = 4;    phi4p = 5; s4p = 6; bt4p = 7;
BC = 3; CD = 2;

%q and its time derivatives at a given t
%q0p=q, q1p = dq0p/dt, q2p = dq1p/dt ...
q0p = [th0p; phi0p; s0p; bt0p];
q1p = [th1p; phi1p; s1p; bt1p];
q2p = [th2p; phi2p; s2p; bt2p];
q3p = [th3p; phi3p; s3p; bt3p];
q4p = [th4p; phi4p; s4p; bt4p];

parameters = [BC,CD];

%BFF method when a basis is given (see function 'basisCD' below) 
[w,a,jk,js] = angularKinQ14(@basisCD,q4p,q3p,q2p,q1p,q0p,parameters);
disp("AKQ's for Example 3 using the BFF method when a basis is given")

disp([w,a,jk,js])

%BFF method when 3 points are given
%In this case, a function for the three points is provided (see function 
%'points123' below). The derivatives are not explicitly provided, but they 
%are calculated using the KinQD04 function.
%
%computing the time derivatives of points123(q,parameters)
[R0p,R1p,R2p,R3p,R4p] = ...
    KinQD04(@points123,q4p,q3p,q2p,q1p,q0p,parameters);

%computing angular kinematic quantities
[w,a,jk,js] = angularKinQ14_3pts(R4p,R3p,R2p,R1p,R0p);
disp("AKQ's for Example 3 using the BFF method when 3 points are given")
disp([w,a,jk,js])

%vector method
%using the same R0p,R1p,R2p,R3p,R4p, values computed above
w = ang_vel_3pts_vec(R1p,R0p);
a = ang_accel_3pts_vec(R2p,R1p,R0p);       
jk = ang_jerk_3pts_vec(R3p,R2p,R1p,R0p);        
js = ang_jounce_snap_3pts_vec(R4p,R3p,R2p,R1p,R0p);
disp("AKQ's for Example 3 using the vector method")
disp([w,a,jk,js])


%--------------------------------------------------------------------------
%auxiliar functions for the above example 
%
%A = points123(q,pars)
%Position vectors of three non-collinear points corresponding to the 
%link CD. These vectors represent rigid extensions of the link and 
%are arranged as columns in matrix A.
%q: generalized coordinates
%pars: parameters independent of time
function fr = points123(q,pars)
  th = q(1); phi = q(2); s=q(3); bt = q(4);

  BC = pars(1); CD = pars(2);

  T1 = HTM(th,[0,0,1],[0,0,0]);
  T2 = HTM(phi,[1,0,0],[s,0,0]);
  T3 = HTM(0,[0,0,1],[0,0,BC]);
  T4 = HTM(bt,[0,0,1],[0,0,0]);
  T5 = HTM(0,[1,0,0],[CD,0,0]);

  RD = T1*T2*T3*T4*T5; RC = T1*T2*T3; RB = T1*T2;
  pD = RD(1:3,4); pC = RC(1:3,4); pB = RB(1:3,4);
  pCD = pD - pC; pBC = pC - pB; 

  x3 = vuni(pCD); z3 = vuni(pBC); y3 = vuni(cross(z3,x3));

  p1 = pC + x3;
  p2 = pC + y3;
  p3 = pC + z3;

  fr = cat(2,p1,p2,p3);
end

%A = basisCD(q,pars)
%Attached basis to the link CD. The vectors are given as columns of 
%matrix A.
%q: generalized coordinates
%pars: parameters independent of time
function fr = basisCD(q,pars)
  th = q(1); phi = q(2); s=q(3); bt = q(4);

  BC = pars(1); CD = pars(2);

  T1 = HTM(th,[0,0,1],[0,0,0]);
  T2 = HTM(phi,[1,0,0],[s,0,0]);
  T3 = HTM(0,[0,0,1],[0,0,BC]);
  T4 = HTM(bt,[0,0,1],[0,0,0]);
  T5 = HTM(0,[1,0,0],[CD,0,0]);

  RD = T1*T2*T3*T4*T5; RC = T1*T2*T3; RB = T1*T2;
  pD = RD(1:3,4); pC = RC(1:3,4); pB = RB(1:3,4);
  pCD = pD - pC; pBC = pC - pB; 

  x3 = vuni(pCD); z3 = vuni(pBC); y3 = vuni(cross(z3,x3));

  fr = cat(2,x3,y3,z3);
end



