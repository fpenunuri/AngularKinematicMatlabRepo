%[p,v,a,j,js] = KinQD04(f,g,t,pars)
%output variables
%p:  f(g(t),pars) or (d^0/dt^0)f = f
%v:  dp/dt or (d^1/dt^1)f
%a:  dv/dt or (d^2/dt^2)f
%j:  da/dt or (d^3/dt^3)f
%js: dj/dt or (d^4/dt^4)f
%input variables
%f: the function to differentiate. F = f(q, pars) is defined by the user
%and F could be a quantity with any number of indices. Nevertheless it is
%ussually a vector or a matrix. The vector q stores the generalized
%coordinates, and pars is an optional parameter independent of time.
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
%[p,v,a,j,js] = KinQD04(f,q4p,q3p,q2p,q1p,q0p,pars)
%output variables
%p:  f(g(t),pars) or (d^0/dt^0)f = f
%v:  dp/dt or (d^1/dt^1)f
%a:  dv/dt or (d^2/dt^2)f
%j:  da/dt or (d^3/dt^3)f
%js: dj/dt or (d^4/dt^4)f
%input variables
%f: the function to differentiate. F = f(q, pars) is defined by the user
%and F could be a quantity with any number of indices. Nevertheless it is
%ussually a vector or a matrix. The vector q stores the generalized
%coordinates, and pars is an optional parameter independent of time.
%q4p,q3p,q2p,q1p,q0p: vectors that store, the time derivatives of 
%the generalized coordinates. Deppending of how is the f function coded, 
%all the vectors can be row or column vectors.
%qkp=(d/dt)^k q; k = 0,1,2,3,4.
%pars: optional vector storing other parameters independent of time, 
%which is the same argument that appears in the function f

function [pos,vel,accel,jerk,jounce_snap] = KinQD04(varargin)
  if isa(varargin{2}, 'function_handle')
    [pos,vel,accel,jerk,jounce_snap] = KinQD04auxF(varargin{:});
  else
    [pos,vel,accel,jerk,jounce_snap] = KinQD04auxP(varargin{:});
  end 
end

%f: position vector (the function to differentiate with respect t)
%g: vector function storing the generalized coordinates
%g = g([q1(t),q2(t),....qm(t)])
%if f depends explicitly on time, use qm=t
%pars: Other (if any) time-independent parameters"; pars is optional
function [pos,vel,accel,jerk,jounce_snap] = KinQD04auxF(f,g,t0,pars)
  if (~exist('pars','var'))
    taux = dual4(t0,1,0,0,0);
    Xt = g(taux);
    frd4 = f(Xt);
    pos = frd4.f0;
    vel = frd4.f1; 
    accel = frd4.f2;
    jerk  = frd4.f3;
    jounce_snap = frd4.f4;
  else    
    taux = dual4(t0,1,0,0,0);
    Xt = g(taux);
    frd4 = f(Xt,pars);
    pos = frd4.f0;
    vel = frd4.f1; 
    accel = frd4.f2;
    jerk  = frd4.f3;
    jounce_snap = frd4.f4;
  end
end

%f: position vector (the function to differentiate with respect t)
%q4p: d4q/dt4 (fourth derivative)
%q3p: d3q/dt3 (third derivative)
%q2p: d2q/dt2 (second derivative)
%q1p: dq/dt (first derivative)
%q1p = (d/dt)[q1(t0),q2(t0),...qm(t0)]
%q0p: d^0 q/dt^0 = q (0th derivative); 
%q0p = [q1(t0),q2(t0),...qm(t0)] generalized coordinates
%if f depends explicitly on time, use qm=t
%pars: Other (if any) time-independent parameters"; pars is optional
function [pos,vel,accel,jerk,jounce_snap] = ...
      KinQD04auxP(f,q4p,q3p,q2p,q1p,q0p,pars)

  if (~exist('pars','var'))
    taux = dual4(0,1,0,0,0);
    Xt = q0p + taux.*q1p + (taux.^2) .* q2p/2 + (taux.^3) .* q3p/6 + ...
         (taux.^4) .* q4p/24;

    frd4 = f(Xt);
    pos = frd4.f0;
    vel = frd4.f1;
    accel = frd4.f2;
    jerk  = frd4.f3;
    jounce_snap = frd4.f4;
  else    
    taux = dual4(0,1,0,0,0);
    Xt = q0p + taux.*q1p + (taux.^2) .* q2p/2 + (taux.^3) .* q3p/6 + ...
         (taux.^4) .* q4p/24;

    frd4 = f(Xt,pars);
    pos = frd4.f0;
    vel = frd4.f1;
    accel = frd4.f2;
    jerk  = frd4.f3;
    jounce_snap = frd4.f4;
  end
end


