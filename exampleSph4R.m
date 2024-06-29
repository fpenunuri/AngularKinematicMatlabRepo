addpath('./01DualNumbersF;./02AuxiliarDualFun;./03KinematicQuantities')
addpath('./04MechFunctions')

%EXAMPLE SPHERICAL 4R MECHANISM
eta1 = pi/2; phi1 = 0; psi = 0; alpha1 = pi/2; alpha2 = pi/6; 
alpha3 = pi/2; alpha4 = pi/2; beta = 0; gm = 0;

%vd: design variable vector 
vd = [eta1,phi1,psi,alpha1,alpha2,alpha3,alpha4,beta,gm];

q0p = 3*pi/4; %input angle th

q1p = 1; %constant angular velocity; thus, q2p, q3p, q4p are zero
q2p = 0; 
q3p = 0;
q4p = 0;

disp("Angular kinematic quantities (AKQ's) from velocity")
disp('up to the jounce/snap as columns')
%BFF method when a basis is given (see function 'BCoupler' and 
%function 'BOutputL' below) 
[w,a,jk,js] = angularKinQ14(@BCoupler,q4p,q3p,q2p,q1p,q0p,vd);
disp("AKQ's for Example 4 using the BFF method when basis is given")
disp('coupler link case')
disp([w,a,jk,js])

[w,a,jk,js] = angularKinQ14(@BOutputL,q4p,q3p,q2p,q1p,q0p,vd);
disp('output link case')
disp([w,a,jk,js])

disp('--------')
disp('3 points are given')
%disp('3 points are given')
%COUPLER 
%BFF method when 3 points are given
%In this case, a function for the three points is provided (see function 
%'Coupler3p' below). The derivatives are not explicitly provided, but they 
%are calculated using the KinQD04 function.
%
%computing the time derivatives of Coupler3p(q,parameters)
[R0p,R1p,R2p,R3p,R4p] = ...
    KinQD04(@Coupler3p,q4p,q3p,q2p,q1p,q0p,vd);

%computing angular kinematic quantities
[w,a,jk,js] = angularKinQ14_3pts(R4p,R3p,R2p,R1p,R0p);
disp("AKQ's for Example 4 using the BFF method")
disp('coupler link case')
disp([w,a,jk,js])

%vector method
%using the same R0p,R1p,R2p,R3p,R4p, values computed above
w = ang_vel_3pts_vec(R1p,R0p);
a = ang_accel_3pts_vec(R2p,R1p,R0p);       
jk = ang_jerk_3pts_vec(R3p,R2p,R1p,R0p);        
js = ang_jounce_snap_3pts_vec(R4p,R3p,R2p,R1p,R0p);
disp("AKQ's for Example 4 using the vector method")
disp('coupler link case')
disp([w,a,jk,js])

%OUTPUT LINK
%BFF method when 3 points are given
%In this case, a function for the three points is provided (see function 
%'OutputL3p' below). The derivatives are not explicitly provided, but they 
%are calculated using the KinQD04 function.
%
%computing the time derivatives of OutputL3p(q,parameters)
[R0p,R1p,R2p,R3p,R4p] = ...
    KinQD04(@OutputL3p,q4p,q3p,q2p,q1p,q0p,vd);

%computing the angular kinematic quantities
[w,a,jk,js] = angularKinQ14_3pts(R4p,R3p,R2p,R1p,R0p);
disp("AKQ's for Example 4 using the BFF method")
disp('output link case')
disp([w,a,jk,js])

%vector method
%using the same R0p,R1p,R2p,R3p,R4p, values computed above
w = ang_vel_3pts_vec(R1p,R0p);
a = ang_accel_3pts_vec(R2p,R1p,R0p);       
jk = ang_jerk_3pts_vec(R3p,R2p,R1p,R0p);        
js = ang_jounce_snap_3pts_vec(R4p,R3p,R2p,R1p,R0p);
disp("AKQ's for Example 4 using the vector method")
disp('output link case')
disp([w,a,jk,js])

%--------------------------------------------------------------------------
%Auxiliar functions
%Some quantities related to the spherical 4R mechanism are calculated by
%[x1,r2,r3,rg,x4,T23,T43,th4] = sph4r_vars(q,par);
%type 'help sph4r_vars' for more information.
%
%A = BCoupler(q,par)
%Attached basis to the coupler link in the spherical 4R mechanism. The 
%vectors are given as columns of matrix A.
%q: generalized coordinates
%pars: parameters independent of time
function fr = BCoupler(q,par)
  [~,r2,~,~,~,T23,~,~] = sph4r_vars(q,par);
  e3c = vuni(r2);
  e1c = T23;
  e2c = cross(e3c,e1c); %already a unit vector since e3c.e1c=0
  fr = cat(2,e1c,e2c,e3c);
end

%A = BOutputL(q,par)
%Attached basis to the output link in the spherical 4R mechanism. The 
%vectors are given as columns of matrix A.
%q: generalized coordinates
%pars: parameters independent of time
function fr = BOutputL(q,par)
  [~,~,~,~,x4,~,T43,~] = sph4r_vars(q,par);
  e3o = vuni(x4);
  e1o = T43;
  e2o = cross(e3o,e1o);
  fr = cat(2,e1o,e2o,e3o);
end

%A = Coupler3p(q,par)
%Position vectors of three non-collinear points corresponding to the 
%coupler link. These vectors represent rigid extensions of the link and 
%are arranged as columns in matrix A.
%q: generalized coordinates
%pars: parameters independent of time
function fr = Coupler3p(q,par)
  [~,r2,~,~,~,T23,~,~] = sph4r_vars(q,par);
  e3c = vuni(r2);
  e1c = T23;
  e2c = cross(e3c,e1c);
  
  p1 = r2 + e1c;
  p2 = r2 + e2c;
  p3 = r2 + e3c;

  fr = cat(2,p1,p2,p3);
end

%A = OutputL3p(q,par)
%Position vectors of three non-collinear points corresponding to the 
%output link. These vectors represent rigid extensions of the link and 
%are arranged as columns in matrix A.
%q: generalized coordinates
%pars: parameters independent of time
function fr = OutputL3p(q,par)
  [~,~,~,~,x4,~,T43,~] = sph4r_vars(q,par);
  e3o = vuni(x4);
  e1o = T43;
  e2o = cross(e3o,e1o);
  
  p1 = x4 + e1o;
  p2 = x4 + e2o;
  p3 = x4 + e3o;

  fr = cat(2,p1,p2,p3);
end
