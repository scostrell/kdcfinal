%HW 5 Prob1 KDC - following notation and calculations from pages 168-175 in MLS
% 1.1 What are the components of the resulting manipulator inertia matrix M(?) 
%(provide the result for element M11)?
% 1.2 What are the components of the resulting Coriolis matrix C(?,? ?) as 
% defined in the equation (4.23) of the book (provide the result for element C21)?
% 1.3 What are the components of the vector N(?,? ?) (provide the result for element N3 )?

syms l0 l1 l2 r0 r1 r2 m1 m2 m3 real positive
syms th1 th2 th3 Ix1 Iy1 Iz1 Ix2 Iy2 Iz2 Ix3 Iy3 Iz3 dth1 dth2 dth3 ddth1...
    ddth2 ddth3 J1 J2 J3 g real
syms m11 m12 m13 m21 m22 m23 m31 m32 m33 %for manipulator matrix
% syms J1(th1, th2, th3) J2(th1, th2, th3) J3(th1, th2, th3)
%%
% M = sym('M',[3 3]);

% These are given in ex 4.3, p. 172
g0sl1 = [eye(3), [0, 0, r0]';
    0, 0, 0, 1];
g0sl2 = [eye(3), [0, r1, l0]';
    0, 0, 0, 1];
g0sl3 = [eye(3), [0, (l1 + r2), l0]';
    0, 0, 0, 1];

% Twists from ex 4.3
tw1 = [0 0 0 0 0 1]';
tw2 = [0 -l0 0 -1 0 0]';
tw3 = [0 -l0 l1 -1 0 0]';

% Matrix exponentials
e1 = twistExp(tw1,th1);
e2 = twistExp(tw2,th2);
e3 = twistExp(tw3,th3);

% Frame configurations via p. 168
gth_sl1 = e1*g0sl1;
gth_sl2 = e1*e2*g0sl2;
gth_sl3 = e1*e2*e3*g0sl3;

% Find elements of body Jacobians via p. 168
tw1_dag_sl1 = inv(adj(e1*g0sl1))*tw1;

tw1_dag_sl2 = inv(adj(e1*e2*g0sl2))*tw1;
tw2_dag_sl2 = inv(adj(e2*g0sl2))*tw2;

tw1_dag_sl3 = inv(adj(e1*e2*e3*g0sl3))*tw1;
tw2_dag_sl3 = inv(adj(e2*e3*g0sl3))*tw2;
tw3_dag_sl3 = inv(adj(e3*g0sl3))*tw3;

% Calculate Jacobians
J1(th1, th2, th3) = [tw1_dag_sl1, zeros(6,2)];
J2(th1, th2, th3) = [tw1_dag_sl2, tw2_dag_sl2, zeros(6,1)];
J3(th1, th2, th3) = [tw1_dag_sl3, tw2_dag_sl3, tw3_dag_sl3];

% Manipulator matrix via manip.m
M =...
manip(m1, m2, m3, J1, J2, J3, Ix1, Iy1, Iz1, Ix2, Iy2, Iz2, Ix3, Iy3, Iz3);
Mnew = M(th1, th2, th3); % to get around indexing issues, see https://www.mathworks.com/matlabcentral/answers/410630-how-to-call-element-of-matrix-of-symbolic-variables
% Coriolis matrix via coriolis.m
C = coriolis(Mnew, th1, th2, th3, dth1, dth2, dth3);
% N via dV.m
N = dV(th1, th2, th3, m1, m2, m3, gth_sl1, gth_sl2, gth_sl3, g);
%% Check with book output
disp("If this is 0, M11 is correct:")
simplify(Iy2*sin(th2)^2 + Iy3*sin(th2+th3)^2 + Iz1 + Iz2*cos(th2)^2+...
    Iz3*cos(th2+th3)^2 + m2*r1^2*cos(th2)^2 + m3*(l1*cos(th2) + ...
    r2*cos(th2+th3))^2 - Mnew(1,1))
%%
disp("If this is 0, C21 is correct:")
simplify(((Iz2 - Iy2 + m2*r1^2)*cos(th2)*sin(th2) + (Iz3 - Iy3)*cos(th2+th3)...
    *sin(th2+th3) + m3*(l1*cos(th2) + r2*cos(th2 + th3))*(l1*sin(th2) + ...
    r2*sin(th2 + th3)))*dth1-C(2,1))
 %%   
disp("If this is 0, N3 is correct:") 
simplify(-m3*g*r2*cos(th2 + th3) - N(3))
%%
disp("M11 is:")
simplify(Mnew(1,1))
disp("C21 is:")
simplify(C(2,1))
disp("N3 is:")
simplify(N(3))

%% Locations of joints

v1 = 2*gth_sl1(1:3,4);
v2 = v1 + 2*(gth_sl2(1:3,4) - v1);
v3 = v2 + 2*(gth_sl3(1:3,4) - v2);

xp(l0, l1, l2, r0, r1, r2, th1, th2, th3) = [0, v1(1), v2(1), v3(1)];
yp(l0, l1, l2, r0, r1, r2, th1, th2, th3) = [0, v1(2), v2(2), v3(2)];
zp(l0, l1, l2, r0, r1, r2, th1, th2, th3) = [0, v1(3), v2(3), v3(3)];
 %% Set values for simulation dynamics
l0_s = 1;
l1_s = 1;
l2_s = 1;
r0_s = l0_s*0.5;
r1_s = l1_s*0.5;
r2_s = l2_s*0.5;
th1_s = 0;
th2_s = 0;
th3_s = 0;
m1_s = 1;
m2_s = 1;
m3_s = 1;
% Via https://en.wikipedia.org/wiki/List_of_moments_of_inertia#List_of_3D_inertia_tensors
Ix1_s = 1/12*m1_s*l0_s^2;
Iy1_s = 1/12*m1_s*l0_s^2;
Iz1_s = 0;
Ix2_s = 1/12*m2_s*l1_s^2;
Iy2_s = 0;
Iz2_s = 1/12*m2_s*l1_s^2;
Ix3_s = 1/12*m3_s*l2_s^2; 
Iy3_s = 0;
Iz3_s = 1/12*m3_s*l2_s^2; 
g_s = 9.81;
%%
% Initial pose
xq = xp(1,1,1,0.5,0.5,0.5,0,0,0);
yq = yp(1,1,1,0.5,0.5,0.5,0,0,0);
zq = zp(1,1,1,0.5,0.5,0.5,0,0,0);

% Desired pose, with th1=0=th3, th2=-pi/2
xdes = xp(1,1,1,0.5,0.5,0.5,0, -pi/2,0);
ydes = yp(1,1,1,0.5,0.5,0.5,0,-pi/2,0);
zdes = zp(1,1,1,0.5,0.5,0.5,0,-pi/2,0);

% plot3(xdes,ydes,zdes,'-o')

%%
% EOMs with constants for simulation
% tau(l0, l1, l2, r0, r1, r2, m1, m2, m3, Ix1, Iy1, Iz1, Ix2, Iy2, Iz2, Ix3, Iy3, Iz3)...
tau = Mnew*[ddth1; ddth2; ddth3] + C*[dth1; dth2; dth3] + N;
tau_s = subs(tau,[l0, l1, l2, r0, r1, r2, m1, m2, m3, Ix1, Iy1, Iz1, Ix2, Iy2, Iz2, Ix3, Iy3, Iz3, g],[l0_s, l1_s, l2_s, r0_s, r1_s, r2_s, m1_s, m2_s, m3_s, Ix1_s,...
    Iy1_s, Iz1_s, Ix2_s, Iy2_s, Iz2_s, Ix3_s, Iy3_s, Iz3_s, g_s]);
tau_update(th1, th2, th3, dth1, dth2, dth3, ddth1, ddth2, ddth3) = tau_s;

M_s = subs(Mnew,[l0, l1, l2, r0, r1, r2, m1, m2, m3, Ix1, Iy1, Iz1, Ix2, Iy2, Iz2, Ix3, Iy3, Iz3],[l0_s, l1_s, l2_s, r0_s, r1_s, r2_s, m1_s, m2_s, m3_s, Ix1_s,...
    Iy1_s, Iz1_s, Ix2_s, Iy2_s, Iz2_s, Ix3_s, Iy3_s, Iz3_s]);
M_update(th1, th2, th3) = M_s;

C_s = subs(C,[l0, l1, l2, r0, r1, r2, m1, m2, m3, Ix1, Iy1, Iz1, Ix2, Iy2, Iz2, Ix3, Iy3, Iz3, g],[l0_s, l1_s, l2_s, r0_s, r1_s, r2_s, m1_s, m2_s, m3_s, Ix1_s,...
    Iy1_s, Iz1_s, Ix2_s, Iy2_s, Iz2_s, Ix3_s, Iy3_s, Iz3_s, g_s]);
C_update(th1, th2, th3, dth1, dth2, dth3) = C_s;

N_s = subs(N,[l0, l1, l2, r0, r1, r2, m1, m2, m3,g],[l0_s, l1_s, l2_s, r0_s, r1_s, r2_s, m1_s, m2_s, m3_s, g_s]);
N_update(th1, th2, th3) = N_s;
% tau_dth(dth1, dth2, dth3, ddth1, ddth2, ddth3) = tau_new;

%% Control
Kp = 106*eye(3);
Kv = 0.001*eye(3);
% Kv = 0.001*eye(3);
% Ki = 0.2*eye(3);
dt = 0.05;
theta_prev = [0, 0, 0]';
dtheta_prev = [0, 0, 0]';
% ddtheta_prev = [0, 0, 0]';
theta_des = [0, -pi/10, 0]';
tau = [0 -0.1 0]';
theta_arr = [];
tau_arr = [];
t = 100;
tol = 0.12;

for i = 1:t
    if (norm(theta_curr - theta_des) < tol) & (theta_des(2) > -pi/2)
        theta_des = theta_des + [0, -pi/10, 0]';
    end
    % Euler integration
    theta_arr = [theta_arr; theta_prev'];
    tau_arr = [tau_arr; tau'];
    ddtheta_curr = double(inv(M_update(theta_prev(1), theta_prev(2), theta_prev(3)))*...
        (tau - N_update(theta_prev(1), theta_prev(2), theta_prev(3))...
        - C_update(theta_prev(1), theta_prev(2), theta_prev(3), dtheta_prev(1),...
        dtheta_prev(2), dtheta_prev(3))*dtheta_prev));
    dtheta_curr = dtheta_prev + dt*ddtheta_curr;
    theta_curr = theta_prev + dt*dtheta_curr;
    % Update error
    err_prev = theta_prev - theta_des;
    err_curr = theta_curr - theta_des;
    d_err = (err_curr - err_prev)/dt;
    tau = -Kv*d_err - Kp*err_curr;
    % Update state
    theta_prev = theta_curr;
    dtheta_prev = dtheta_curr;
    ddtheta_prev = ddtheta_curr;
end

figure
hold on
plot(1:t,theta_arr(:,1))
plot(1:t,theta_arr(:,2))
plot(1:t,theta_arr(:,3))
legend('\theta_1','\theta_2','\theta_3')
hold off

figure
hold on
plot(1:t,tau_arr(:,1))
plot(1:t,tau_arr(:,2))
plot(1:t,tau_arr(:,3))
legend('\tau_1','\tau_2','\tau_3')
hold off

% theta = [th1, th2, th3]';

% theta_prev = [0, 0, 0]';
% theta_curr = [0.1, 0.1, 0.1]';
% dtheta = (theta_curr - theta_prev)/dt;
% ddtheta = [0,0,0]';
% err_prev = theta_prev - theta_des;
% err_curr = theta_curr - theta_des;
% d_err = (err_curr - err_prev)/dt;

% tau_in = -Kv*d_err - Kp*err_curr;
% tau_s = tau_dth(dtheta(1), dtheta(2), dtheta(3), ddtheta(1), ddtheta(2), ddtheta(3));

% dtheta = (theta_des - theta_curr);
% ddtheta = [0,0.1,0]';
% tau_des = tau_fullsolve(theta_des(1),theta_des(2),theta_des(3),dtheta(1), dtheta(2), dtheta(3), ddtheta(1), ddtheta(2), ddtheta(3));
% solve_for_th = tau_dth(dtheta(1), dtheta(2), dtheta(3), ddtheta(1), ddtheta(2), ddtheta(3));

%%
% eqns = [tau_des(1) == solve_for_th(1), tau_des(2) == solve_for_th(2), tau_des(3) == solve_for_th(3)];
% S = solve(eqns,[th1,th2,th3]);


