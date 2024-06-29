%w = ang_vel_3pts_vec(R1p,R0p)
%output variables
%w: angular velocity
%input variables
%R1p,R0p: matrices that store, as columns, the velocity and position, 
%respectively, of three non-collinear points.

function fr = ang_vel_3pts_vec(RP,R)
  r1 = R(:,1); r2 = R(:,2);
  v1 = RP(:,1); v2 = RP(:,2); v3 = RP(:,3);

  a = v3 - v1; b = v2 - v1; c = r2 - r1;
  ac = sum(a.*c);

  if all(abs(ac)<=eps)
    error('division by 0, use ang_vel_3pts instead')
  end

  fr = Scross(a)*b./ac;
end
