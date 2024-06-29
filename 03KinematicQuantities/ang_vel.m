%w = ang_vel(fBase,g,t,pars)
%output variables
%w: angular velocity
%input variables
%fBase: F = fBase(q, pars) returns a matrix F, where each column represents 
%a basis vector attached to the rigid body. This function depends on the 
%generalized coordinates q and the optional parameter pars, which remains 
%independent of time.
%g: q = g(t) returns [q1(t);q2(t);...;qm(t)] or [q1(t),q2(t),...,qm(t)], 
%the vector of generalized coordinates as function of t. If f deppends 
%explicitly on time, use qm(t)=t. The position vector as function of t is 
%r = f(g(t),pars)
%t: the time of evaluation
%pars: optional vector storing other parameters independent of time, 
%which is the same argument that appears in the function fBase.
%This function utilizes the BFF method and is particularly useful for 
%generating a "continuous" trajectory (r: R --> R^3) for the angular 
%velocity.
%
%Alternatively, this function can also be invoked as:
%w = ang_vel(fBase,q1p,q0p,pars)
%output variables
%w: angular velocity
%input variables
%fBase: F = fBase(q, pars) returns a matrix F, where each column represents 
%a basis vector attached to the rigid body. This function depends on the 
%generalized coordinates q and the optional parameter pars, which remains 
%independent of time.
%q1p,q0p: column matrices that store, the time derivatives of the 
%generalized coordinates
%qkp=(d/dt)^k q; k = 0,1.
%pars: optional vector storing other parameters independent of time, 
%which is the same argument that appears in the function fBase.
%This function uses the BFF method.

function fr = ang_vel(varargin)
  [pos,vel] = KinQD01(varargin{:});    
  e1 = pos(:,1); e2 = pos(:,2); e3 = pos(:,3);
  e1v = vel(:,1); e2v = vel(:,2); e3v = vel(:,3);    
  fr = sp(e2v,e3)*e1 + sp(e3v,e1)*e2 + sp(e1v,e2)*e3;
end 
