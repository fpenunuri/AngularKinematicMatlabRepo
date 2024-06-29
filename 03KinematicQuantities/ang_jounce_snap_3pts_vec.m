%ajs = ang_jounce_snap_3pts_vec(R4p,R3p,R2p,R1p,R0p)
%output variables
%ajs: angular jounce/snap
%input variables
%R4p,R3p,R2p,R1p,R0p: matrices that store, as columns, the jounce/snap, 
%jerk, acceleration, velocity and position, respectively, of three 
%non-collinear points.

function fr = ang_jounce_snap_3pts_vec(R4p,R3p,R2p,R1p,R0p)

  r1  = R0p(:,1); r2  = R0p(:,2);
  v1  = R1p(:,1); v2  = R1p(:,2); v3  = R1p(:,3);
  a1  = R2p(:,1); a2  = R2p(:,2); a3  = R2p(:,3);
  jk1 = R3p(:,1); jk2 = R3p(:,2); jk3 = R3p(:,3);
  js1 = R4p(:,1); js2 = R4p(:,2); js3 = R4p(:,3);

  a   = v3  - v1;  b   = v2  - v1;  c   = r2  - r1;
  a1p = a3  - a1;  b1p = a2  - a1;  c1p = v2  - v1;
  a2p = jk3 - jk1; b2p = jk2 - jk1; c2p = a2  - a1;
  a3p = js3 - js1; b3p = js2 - js1; c3p = jk2 - jk1;

  ac = sp(a,c);
  if all(abs(ac)<=eps)
    error('division by 0, use angularKinQ14_3pts instead')
  end

  apcacp = sp(a1p,c) + sp(a,c1p);

  num1 = Scross(a3p)*b + 3*Scross(a2p)*b1p + 3*Scross(a1p)*b2p + ...
         Scross(a)*b3p;

  num2 = 3*(Scross(a2p)*b+2*Scross(a1p)*b1p + Scross(a)*b2p)*apcacp;

  num3 = 3*(Scross(a1p)*b + Scross(a)*b1p)*(sp(a2p,c) + 2*sp(a1p,c1p) + ...
                                            sp(a,c2p));

  num4 = 6*(Scross(a1p)*b + Scross(a)*b1p)*apcacp.^2;

  num5 = Scross(a)*b .*(sp(a3p,c) + 3*sp(a2p,c1p) + 3*sp(a1p,c2p) + ...
                        sp(a,c3p));

  num6 = 6*Scross(a)*b .*(sp(a2p,c) + 2*sp(a1p,c1p) + sp(a,c2p)).*apcacp; 

  num7 = 6*Scross(a)*b.*apcacp.^3;

  fr = num1./ac - num2./ac.^2 - num3./ac.^2 + num4./ac.^3 - num5./ac.^2 + ...
       num6./ac.^3 - num7./ac.^4;
end
