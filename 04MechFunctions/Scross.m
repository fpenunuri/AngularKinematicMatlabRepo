%Matrix representation (in the standard basis of R3) for the cross product 
%operator
%a x v = Scross(a) * v
function fr = Scross(a)
  fr = [0,-a(3),a(2);a(3),0,-a(1);-a(2),a(1),0];
end


