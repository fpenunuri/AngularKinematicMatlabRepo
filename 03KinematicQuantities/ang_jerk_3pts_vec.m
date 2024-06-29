%aj = ang_jerk_3pts_vec(R3p,R2P,R1P,R0p)
%output variables
%aj: angular jerk
%input variables
%R3p,R2p,R1p,R0p: matrices that store, as columns, the jerk, acceleration, 
%velocity and position, respectively, of three non-collinear points.

function fr = ang_jerk_3pts_vec(R3p,R2p,R1p,R0p)
  r1  = R0p(:,1); r2  = R0p(:,2);
  v1  = R1p(:,1); v2  = R1p(:,2); v3  = R1p(:,3);
  a1  = R2p(:,1); a2  = R2p(:,2); a3  = R2p(:,3);
  jk1 = R3p(:,1); jk2 = R3p(:,2); jk3 = R3p(:,3);

  a   = v3  - v1;  b   = v2  - v1;  c   = r2 - r1;
  a1p = a3  - a1;  b1p = a2  - a1;  c1p = v2 - v1;
  a2p = jk3 - jk1; b2p = jk2 - jk1; c2p = a2 - a1;

  ac = sp(a,c);

  if all(abs(ac)<=eps)
    error('division by 0, use angularKinQ13_3pts instead')
  end

  term1 = (Scross(a2p)*b + 2*Scross(a1p)*b1p + Scross(a)*b2p);

  term2 = 2*(Scross(a1p)*b + Scross(a)*b1p)*(sp(a1p,c) + sp(a,c1p)) + ...
          Scross(a)*b .*(sp(a2p,c) + 2*sp(a1p,c1p) + sp(a,c2p));

  term3 =  2* Scross(a)*b .* (sp(a1p,c) + sp(a,c1p)).^2;

  fr = term1./ac - term2./(ac.^2) + term3./(ac.^3);
end
