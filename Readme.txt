The functions in folder "01DualNumbersF" were downloaded from:
R. Peón-Escalante, K.B. Cantún-Avila, O. Carvente, A. Espinosa-Romero, F.
Peñuñuri, Dual number implementation to compute higher order directional
derivatives, Mendeley Data, v1, 2022. http://dx.doi.org/10.17632/kcrm6pmk7d.
3.

We have added the "zerosLike" function to the "dualk.m" files (k=1,2,3,4).

It is also recommended ref:
R. Peón-Escalante, K. Cantún-Avila, O. Carvente, A. Espinosa-Romero,
F. Peñuñuri, A dual number formulation to efficiently compute higher order
directional derivatives, Journal of Computational Science 76 (2024) 102217.
doi:10.1016/j.jocs.2024.102217.


The folders "02AuxiliarDualFun", "03KinematicQuantities" and "04MechFunctions"
contain the implementation of our study about angular kinematic quantities.

The files "RCR_a.m", "RCR_b.m", "exampleSph4R.m" and "inv_acc_3pts" contain the
examples of the paper.

For help, on a matlab command windows type: help function_name, where
function_name is one of: 
a) folder "02AuxiliarDualFun":
args2duals, f2dualf

b) folder "03KinematicQuantities":
KinQD01, KinQD02, KinQD03, KinQD04, ang_vel, angularKinQ12, angularKinQ13,
angularKinQ14, ang_vel_3pts, angularKinQ12_3pts, angularKinQ13_3pts,
angularKinQ14_3pts, ang_vel_3pts_vec, angularKinQ12_3pts_vec,
angularKinQ13_3pts_vec, angularKinQ14_3pts_vec

c) folder "04MechFunctions":
absX, HTM, rot_mat, Scross, sp, sph4r_vars, vtangent, vuni


Although care has been taken on writing the Matlab code, 
there is not guaranteed to be free of error. The use of the dual 
classes it is entirely at the user's risk and the author disclaims all 
liability.

F. Peñuñuri
