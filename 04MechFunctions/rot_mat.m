%fr = rot_mat(th,n) is the rotation matrix of angle th around axis n
%n is not necessarily a unit vector
function fr = rot_mat(th,eje)     
  %making 'eje' a unit vector
  ejeu = eje./sqrt(sum(eje.*eje));
  
  n1 = ejeu(1); n2 = ejeu(2); n3 = ejeu(3);
  
  fr(1,:) = [1 + (n2.^2 + n3.^2).*(cos(th) - 1), n1.*n2 - n1.*n2.* ...
            cos(th) - n3.*sin(th), n1.*n3 - n1.*n3.*cos(th) + n2.*sin(th)];

  fr(2,:) = [n1.*n2 - n1.*n2.*cos(th) + n3.*sin(th), 1 + (n1.^2 + ...
            n3.^2).*(cos(th) - 1), n2.*n3 - n2.*n3.*cos(th) - n1.*sin(th)];

  fr(3,:) = [n1.*n3 - n1.*n3.*cos(th) - n2.*sin(th), n2.*n3 - n2.*n3.* ...
             cos(th) + n1.*sin(th), 1 + (n1.^2 + n2.^2).*(cos(th) - 1)];
end
