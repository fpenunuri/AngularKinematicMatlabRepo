%[p,v,a,j] = KinQD03(f,g,t,pars)
%output variables
%p:  f(g(t),pars) or (d^0/dt^0)f = f
%v:  dp/dt or (d^1/dt^1)f
%a:  dv/dt or (d^2/dt^2)f
%j:  da/dt or (d^3/dt^3)f
%input variables
%f: the function to differentiate. F = f(q, pars), this function is defined 
%by the user and can represent a quantity with any number of indices. 
%However, it is typically a vector or a matrix. The vector q stores the 
%generalized coordinates, while pars is an optional parameter independent 
%of time.
%g: q = g(t) returns the vector of generalized coordinates as a function 
%of t. Depending on how the f function is defined, q could be a row or a 
%column vector. If f explicitly depends on time, use qm(t) = t. 
%t: the time of evaluation.
%pars: optional vector storing other parameters independent of time, 
%which is the same argument that appears in the function f.
%If f is a position vector, this function is particularly useful for 
%generating a "continuous" trajectory (r: R --> R^3) for the  kinematic 
%quantities.
%
%Alternatively, this function can also be invoked as:
%[p,v,a,j] = KinQD03(f,q3p,q2p,q1p,q0p,pars)
%output variables
%p:  f(g(t),pars) or (d^0/dt^0)f = f
%v:  dp/dt or (d^1/dt^1)f
%a:  dv/dt or (d^2/dt^2)f
%j:  da/dt or (d^3/dt^3)f
%input variables
%f: the function to differentiate. F = f(q, pars), F is defined by the user
%it could be a quantity with any number of indices. Nevertheless ussually 
%is a vector or a matrix. The vector q stores the generalized coordinates, 
%pars is an optional parameter independent of time.
%q3p,q2p,q1p,q0p: vectors that store, the time derivatives of 
%the generalized coordinates. Deppending of how is the f function coded, 
%all the vectors can be row or column vectors.
%qkp=(d/dt)^k q; k = 0,1,2,3.
%pars: optional vector storing other parameters independent of time, 
%which is the same argument that appears in the function f

function [pos,vel,accel,jerk] = KinQD03(varargin)
  if isa(varargin{2}, 'function_handle')
    [pos,vel,accel,jerk] = KinQD03auxF(varargin{:});
  else
    [pos,vel,accel,jerk] = KinQD03auxP(varargin{:});
  end 
end

%f: position vector (the function to differentiate with respect t)
%g: vector function storing the generalized coordinates
%g = g([q1(t),q2(t),....qm(t)])
%if f depends explicitly on time, use qm=t
%pars: Other (if any) time-independent parameters"; pars is optional
function [pos,vel,accel,jerk] = KinQD03auxF(f,g,t0,pars)
  if (~exist('pars','var'))
    taux = dual3(t0,1,0,0);
    Xt = g(taux);
    frd3 = f(Xt);
    pos = frd3.f0;
    vel = frd3.f1; 
    accel = frd3.f2;
    jerk  = frd3.f3;    
  else    
    taux = dual3(t0,1,0,0);
    Xt = g(taux);
    frd3 = f(Xt,pars);
    pos = frd3.f0;
    vel = frd3.f1; 
    accel = frd3.f2;
    jerk  = frd3.f3;
  end
end


%f: position vector (the function to differentiate with respect t)
%q3p: d3q/dt3 (third derivative)
%q2p: d2q/dt2 (second derivative)
%q1p: dq/dt (first derivative)
%q1p = (d/dt)[q1(t0),q2(t0),...qm(t0)]
%q0p: d^0 q/dt^0 = q (0th derivative); 
%q0p = [q1(t0),q2(t0),...qm(t0)] generalized coordinates
%if f depends explicitly on time, use qm=t
%pars: Other (if any) time-independent parameters"; pars is optional
function [pos,vel,accel,jerk] = KinQD03auxP(f,q3p,q2p,q1p,q0p,pars)
  if (~exist('pars','var'))
    taux = dual3(0,1,0,0);
    Xt = q0p + taux.*q1p + (taux.^2) .* q2p/2 + (taux.^3) .* q3p/6;
    frd3 = f(Xt);
    pos = frd3.f0;
    vel = frd3.f1;
    accel = frd3.f2;
    jerk  = frd3.f3;
  else   
    taux = dual3(0,1,0,0);
    Xt = q0p + taux.*q1p + (taux.^2) .* q2p/2 + (taux.^3) .* q3p/6;
    frd3 = f(Xt,pars);
    pos = frd3.f0;
    vel = frd3.f1;
    accel = frd3.f2;
    jerk  = frd3.f3;
  end
end

