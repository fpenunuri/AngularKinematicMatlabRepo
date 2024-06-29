addpath('./01DualNumbersF;./02AuxiliarDualFun;./03KinematicQuantities')
addpath('./04MechFunctions')

%The RCR robot manipulator
%Angular Kinematic quantities, given the generalized coordinates as
%functions of time

BC = 3; CD = 2;
parameters = [BC,CD];

%Example of generalized coordinates as a function of time
fq = @(t) [cos(t),sin(t),sin(t).*cos(t),sin(cos(t))]';


%AKQs as a function of t
AKQs = @(t) angularKinQ14(@basisCD,fq,t,parameters);

%initializing matrices w, a, jk, and js to plot data
np=100;
tvec = linspace(0,2*pi,np);
w = zeros(3,np);
a = zeros(3,np);
jk = zeros(3,np);
js = zeros(3,np);

%storing values in matrices
for k=1:np
    [w(:,k),a(:,k),jk(:,k),js(:,k)]=AKQs(tvec(k));
end

close all;

%new figure
figure;

%subplot 1
subplot(2, 2, 1);
x = w(1,:);
y = w(2,:);
z = w(3,:);
plot3(x, y, z, 'Color', [1,0,0], 'LineWidth', 2);
title('Angular velocity');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;

%subplot 2
subplot(2, 2, 2);
x = a(1,:);
y = a(2,:);
z = a(3,:);
plot3(x, y, z, 'Color', [1,0,1]/sqrt(2), 'LineWidth', 2);
title('Angular acceleration');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;

%subplot 3
subplot(2, 2, 3);
x = jk(1,:);
y = jk(2,:);
z = jk(3,:);
plot3(x, y, z, 'Color', [0.8937,0.4486,0], 'LineWidth', 2);
title('Angular jerk');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;

%subplot 4
subplot(2, 2, 4);
x = js(1,:);
y = js(2,:);
z = js(3,:);
plot3(x, y, z, 'Color', [0,0,1], 'LineWidth', 2);
title('Angular jounce/snap');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;

%A = basisCD(q,pars)
%Attached basis to the link CD. The vectors are given as columns of 
%matrix A.
%q: generalized coordinates
%pars: parameters independent of time
function fr = basisCD(q,pars)
  th = q(1); phi = q(2); s=q(3); bt = q(4);

  BC = pars(1); CD = pars(2);

  T1 = HTM(th,[0,0,1],[0,0,0]);
  T2 = HTM(phi,[1,0,0],[s,0,0]);
  T3 = HTM(0,[0,0,1],[0,0,BC]);
  T4 = HTM(bt,[0,0,1],[0,0,0]);
  T5 = HTM(0,[1,0,0],[CD,0,0]);

  RD = T1*T2*T3*T4*T5; RC = T1*T2*T3; RB = T1*T2;
  pD = RD(1:3,4); pC = RC(1:3,4); pB = RB(1:3,4);
  pCD = pD - pC; pBC = pC - pB; 

  x3 = vuni(pCD); z3 = vuni(pBC); y3 = vuni(cross(z3,x3));

  fr = cat(2,x3,y3,z3);
end
