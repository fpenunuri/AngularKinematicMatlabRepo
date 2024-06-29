%[x1,r2,r3,rg,x4,T23,T43,th4] = sph4r_vars(thv,vd) gives some 
%quantities related to the spherical 4R mechanism
%OUTPUT VARIABLES
%x1: position for the crank-frame joint
%r2: position for the crank-coupler joint
%r3: position for the coupler-oscillator joint
%rg:coupler point vector for the spherical 4R mechanism
%x4: position for the frame-oscillator joint
%T23: unit tangent vector in r2 to the arc from r2 to r3
%T43: unit tangent vector in x4 to the arc from x4 to r3
%th4: output angle
%INPUT VARIABLES
%th: input angle
%vd: design variable vector 
%vd = [eta1,phi1,psi,alpha1,alpha2,alpha3,alpha4,beta,gm] 
%a small quantity is taken as 1e-13 
function  [x1,r2,r3,rg,x4,T23,T43,th4] = sph4r_vars(thv,vd)
  vdc = num2cell(vd,1);
  [eta1,phi1,psi,alpha1,alpha2,alpha3,alpha4,beta,gm] = vdc{:};
  th = thv(1);
  
  epsilon = 1e-13;
  if(eta1 == 0)
    eta_aux = epsilon;
  elseif(eta1 == pi) 
    eta_aux = pi - epsilon;
  else
    eta_aux = eta1;
  end
  
  E3  =eye(3);
  e1 = E3(:,1); e2 = E3(:,2); e3 = E3(:,3);
  
  x1 = cos(phi1).*sin(eta_aux).*e1 + sin(phi1).*sin(eta_aux).*e2 + ...
       cos(eta_aux).*e3;
  
  x1p = x1 - cos(eta_aux).*e3;
  x1T = rot_mat(pi/10.0, e3) * x1p;
  T1 = vtangent(x1p, x1T);
  T14 = rot_mat(psi, x1) * T1;
  n14 = vuni(cross(x1, T14));
  x4 = rot_mat(alpha1, n14) * x1;
  r1 = rot_mat(alpha2, n14) * x1;
  r2 = rot_mat(th, x1) * r1;
  A = sin(alpha2).*sin(alpha4).*sin(th);    
  B = cos(alpha2).*sin(alpha4).*sin(alpha1) - ...
      sin(alpha2).*sin(alpha4).*cos(alpha1).*cos(th);

  C = cos(alpha2).*cos(alpha4).*cos(alpha1) + ...
      sin(alpha2).*cos(alpha4).*sin(alpha1).*cos(th) - ...
      cos(alpha3);

  th4 = 2*atan((-A - sqrt(A.^2 + B.^2 - C.^2))./(C - B));
  r3 = rot_mat(th4, -x4) * rot_mat(alpha4, -n14) * x4;
  n23 = vuni(cross(r2, r3));
  rb = rot_mat(beta, n23) * r2;
  rbg = rot_mat(beta + gm, n23) * r2;

  auxv = absX(acos(sp(r2,r3))-alpha3);
  logicval = auxv <= epsilon;

  if(logicval)
    rg = rot_mat(pi/2.0, rb) * rbg;
    T23 = vtangent(r2, r3);
    T43 = vtangent(x4, r3);
  else
    rg = NaN.*e1; T23 = NaN.*e1; T43 = NaN.*e1;
    th4 = NaN.*e1; r2 = NaN.*e1; r3 = NaN.*e1;
  end
end
