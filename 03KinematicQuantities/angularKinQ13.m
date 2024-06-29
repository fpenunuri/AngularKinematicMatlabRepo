%[w,aa,aj] = angularKinQ13(fBase,g,t,pars)
%output variables
%w: angular velocity
%aa: angular acceleration
%aj: angular jerk
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
%kinematic quantities.
%
%Alternatively, this function can also be invoked as:
%[w,aa,aj] = angularKinQ13(fBase,q3p,q2p,q1p,q0p,pars)
%output variables
%w: angular velocity
%aa: angular acceleration
%aj: angular jerk
%input variables
%fBase: F = fBase(q, pars) returns a matrix F, where each column represents 
%a basis vector attached to the rigid body. This function depends on the 
%generalized coordinates q and the optional parameter pars, which remains 
%independent of time.
%q3p,q2p,q1p,q0p: column matrices that store, the time derivatives of 
%the generalized coordinates
%qkp=(d/dt)^k q; k = 0,1,2,3.
%pars: optional vector storing other parameters independent of time, 
%which is the same argument that appears in the function fBase
%This function uses the BFF method.

function [omega,alpha,ang_jerk] = angularKinQ13(varargin)

  [pos,vel,accel,jerk] = KinQD03(varargin{:});
  
  pc = num2cell(pos,1); vc = num2cell(vel,1); ac = num2cell(accel,1);
  jc = num2cell(jerk,1);
  
  [e1,e2,e3]    = pc{:};
  [e1v,e2v,e3v] = vc{:};
  [e1a,e2a,e3a] = ac{:};
  [e1j,e2j,e3j] = jc{:};
  
  ang_jerk = (2*sp(e2a,e3v) + sp(e2j,e3) + sp(e2v,e3a)).*e1 + ...
      (2*sp(e3a,e1v) + sp(e3j,e1) + sp(e3v,e1a)).*e2 + ...
      (2*sp(e1a,e2v) + sp(e1j,e2) + sp(e1v,e2a)).*e3 + ...
      2*(sp(e2a,e3) + sp(e2v,e3v)).*e1v + sp(e2v,e3).*e1a + ...
      2*(sp(e3a,e1) + sp(e3v,e1v)).*e2v + sp(e3v,e1).*e2a + ...
      2*(sp(e1a,e2) + sp(e1v,e2v)).*e3v + sp(e1v,e2).*e3a;  
  
  alpha = (sp(e2a,e3)+sp(e2v,e3v)).*e1 + sp(e2v,e3).*e1v + ...
          (sp(e3a,e1)+sp(e3v,e1v)).*e2 + sp(e3v,e1).*e2v + ...
          (sp(e1a,e2)+sp(e1v,e2v)).*e3 + sp(e1v,e2).*e3v; 
  
  omega = sp(e2v,e3).*e1 + sp(e3v,e1).*e2 + sp(e1v,e2).*e3;
end  
