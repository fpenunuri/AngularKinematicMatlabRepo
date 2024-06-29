%fr = HTM(th, n, D) gives the homogeneous transformation matrix of angle th 
%around vector n and displacement D.
%This function might require changing its arguments to dual arguments
%depending on whether some of them are dual. Alternatively, you can use
%the f2dualf function.
%
function fr = HTM(th,eje,Dr)
  vars = args2duals(th,eje,Dr);
  [th,eje,Dr] = vars{:};
  fr = zeros(4,'like',th);

  fr(1:3,1:3) = rot_mat(th,eje);
  fr(1:3,4) = Dr;
  fr(4,:) = [0,0,0,1];
end

