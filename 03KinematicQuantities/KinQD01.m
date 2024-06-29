%[p,v] = KinQD01(f,g,t,pars)
%output variables
%p:  f(g(t),pars) or (d^0/dt^0)f = f
%v:  dp/dt or (d^1/dt^1)f
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
%generating a "continuous" trajectory (r: R --> R^3) for the kinematic 
%quantities.
%
%Alternatively, this function can also be invoked as:
%[p,v] = KinQD01(f,q1p,q0p,pars)
%output variables
%p:  f(g(t),pars) or (d^0/dt^0)f = f
%v:  dp/dt or (d^1/dt^1)f
%input variables
%f: the function to differentiate. F = f(q, pars), F is defined by the user
%it could be a quantity with any number of indices. Nevertheless ussually 
%is a vector or a matrix. The vector q stores the generalized coordinates, 
%pars is an optional parameter independent of time.
%q1p,q0p: vectors that store, the time derivatives of the generalized 
%coordinates. Deppending of how is the f function coded, all the vectors 
%can be row or column vectors.
%qkp=(d/dt)^k q; k = 0,1.
%pars: optional vector storing other parameters independent of time, 
%which is the same argument that appears in the function f

function [pos,vel] = KinQD01(varargin)
  if isa(varargin{2}, 'function_handle')
    [pos,vel] = KinQD01auxF(varargin{:});
  else
    [pos,vel] = KinQD01auxP(varargin{:});
  end 
end

%f: position vector (the function to differentiate with respect t)
%f = f(q,pars), pars is optional
%g: vector function storing the generalized coordinates
%g = g([q1(t),q2(t),....qm(t)])
%if f depends explicitly on time, use qm=t
%pars: Other (if any) time-independent parameters"; pars is optional
function [pos,vel] = KinQD01auxF(f,g,t0,pars)
  if (~exist('pars','var'))
    taux = dual1(t0,1);
    Xt = g(taux);
    frd1 = f(Xt);
    pos = frd1.f0;
    vel = frd1.f1;
  else
    taux = dual1(t0,1);
    Xt = g(taux);
    frd1 = f(Xt,pars);
    pos = frd1.f0;
    vel = frd1.f1; 
  end
end

%f: position vector (the function to differentiate with respect t)
%f = f(q,pars), pars is optional
%q1p: dq/dt (first derivative)
%q1p = (d/dt)[q1(t0),q2(t0),...qm(t0)]
%q0p: d^0 q/dt^0 = q (0th derivative); 
%q0p = [q1(t0),q2(t0),...qm(t0)] generalized coordinates
%if f depends explicitly on time, use qm=t
%pars: Other (if any) time-independent parameters"; pars is optional
function [pos,vel] = KinQD01auxP(f,q1p,q0p,pars)
  if (~exist('pars','var'))
    taux = dual1(0,1);
    Xt = q0p + taux.*q1p;
    frd1 = f(Xt);
    pos = frd1.f0;
    vel = frd1.f1; 
  else
    taux = dual1(0,1);
    Xt = q0p + taux.*q1p;
    frd1 = f(Xt,pars);
    pos = frd1.f0;
    vel = frd1.f1; 
  end
end

