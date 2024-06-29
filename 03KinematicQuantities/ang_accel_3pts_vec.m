%aa = ang_accel_3pts_vec(R2P,R1P,R0p)
%output variables
%aa: angular acceleration
%input variables
%R2p,R1p,R0p: matrices that store, as columns, the acceleration, velocity 
%and position, respectively, of three non-collinear points.

function fr = ang_accel_3pts_vec(R2p,R1p,R0p)
  r1 = R0p(:,1); r2 = R0p(:,2);
  v1 = R1p(:,1); v2 = R1p(:,2); v3 = R1p(:,3);
  a1 = R2p(:,1); a2 = R2p(:,2); a3 = R2p(:,3);

  a   = v3 - v1; b   = v2 - v1;  c   = r2 - r1;
  a1p = a3 - a1; b1p = a2 - a1 ; c1p = v2 - v1;

  ac = sp(a,c);
  if all(abs(ac)<=eps)
    error('division by 0, use angularKinQ12_3pts instead')
  end

  fr = (Scross(a1p)*b + Scross(a)*b1p)./ac - ...
       (sp(a1p,c) + sp(a,c1p)).*Scross(a)*b ./ac.^2;
end

