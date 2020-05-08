% Project EOM for Barrett arm, see p.172 of MLS
%%
syms th1 th2 th3 th4 th5 th6 th7 dth1 dth2 dth3 dth4 dth5 dth6 dth7 real

% First, need to calculate inertia matrices Mi:
Q0 = [0.763, 0.64, -0.089; 
    0.062, 0.065, 0.996;
    0.643, -0.765, 0.01];
Q1 = [-0.004, 0.097, 0.995; 
    0.03, -0.995, 0.097;
    0.999, 0.031, 0.001];
Q2 = [0.025, 0.035, 0.999; 
    0.85, 0.525, -0.039;
    -0.525, 0.85, -0.016];
Q3 = [-0.044, -0.99, -0.132;
    0.999, -0.044, -0.005;
    0, -0.132, 0.991];
Q4 = [-0.123, -0.053, 0.991;
    -0.002, -0.998, -0.053;
    0.992, -0.01, 0.123];
Q5 = [0.999, -0.01, -0.01;
    0.01, 0.247, 0.969;
    0, -0.969, 0.248];
Q6 = [0.002, -0.006, 0.999;
    0.979, 0.204, 0;
    -0.205, 0.979, 0.007];
Q7 = [0.913, 0.408, 0;
    -0.408, 0.913, 0;
    0, 0, 0.999];

m0 = 10.17;
m1 = 10.23;
m2 = 4.06;
m3 = 1.7;
m4 = 2.5;
m5 = 0.12;
m6 = 0.5;
m7 = 0.08;

I0 = [0.089, 0.142, 0.187];
I1 = [0.085, 0.107, 0.128];
I2 = [0.013, 0.018, 0.021];
I3 = [0.003, 0.058, 0.058];
I4 = [0.004, 0.016, 0.016];
I5 = [4.5, 5.4, 6.9]*1e-5;
I6 = [2.7, 5.75, 6.8]*1e-4;
I7 = [4.2, 4.21, 8.45]*1e-5;

tw1 = [0.72 -0.61 0 0 0 1]';
tw2 = [0 -1.346 0.72 -1 0 0]';
tw3 = tw1;
tw4 = [0 -1.896 0.765 -1 0 0]';
tw5 = tw1;
tw6 = [0 -2.196 0.72 -1 0 0]';
tw7 = tw1;
% twists = [tw1 tw2 tw3 tw4 tw5 tw6 tw7];

base = [0.75, 0.5, 1]; % location of base frame in meters
base_rot = [0 -1 0;
            1 0 0;
            0 0 1];
% alt_rot= % need to rotate some frames differently?

% COM of link within link frame
r0 = [-0.015, -0.271, -0.145];
r1 = [-0.007, 0.127, 0];
r2 = [-0.002, 0.031, 0.015];
r3 = [-0.042, 0.210, 0];
r4 = [0.003, 0, 0.138];
r5 = [0, 0.007, 0.002];
r6 = [0, -0.024, 0.028];
r7 = [-1.8, 2.6, 35.3]*1e-4;

% link frame locations
L0 = base;
L1 = [0.61, 0.72, 1.346];
L2 = L1;
L3 = L1;
L4 = [0.61,0.765,1.896];
L5 = L1 + [0, 0, 0.85];
L6 = L5;
L7 = L1 + [0, 0, 0.91];

% marker, spatial to tool location

g0sl0 =  [base_rot,(L0 + r0)';
    0, 0, 0, 1];
g0sl1 = [base_rot, (L1 + r1)';
    0, 0, 0, 1];
g0sl2 = [base_rot, (L2 + r2)';
    0, 0, 0, 1];
g0sl3 = [base_rot, (L3 + r3)';
    0, 0, 0, 1];
g0sl4 = [base_rot, (L4 + r4)';
    0, 0, 0, 1];
g0sl5 = [base_rot, (L5 + r5)';
    0, 0, 0, 1];
g0sl6 = [base_rot, (L6 + r6)';
    0, 0, 0, 1];
g0sl7 = [base_rot, (L7 + r7)';
    0, 0, 0, 1];

% Matrix exponentials
fprintf(1, "Matrix exponentials\n");
e1 = twistExp(tw1,th1);
e2 = twistExp(tw2,th2);
e3 = twistExp(tw3,th3);
e4 = twistExp(tw4,th4);
e5 = twistExp(tw5,th5);
e6 = twistExp(tw6,th6);
e7 = twistExp(tw7,th7);

% Frame configurations via MLS p. 168
fprintf(1, "Frame configurations\n");
gth_sl1 = e1*g0sl1;
gth_sl2 = e1*e2*g0sl2;
gth_sl3 = e1*e2*e3*g0sl3;
gth_sl4 = e1*e2*e3*e4*g0sl4;
gth_sl5 = e1*e2*e3*e4*e5*g0sl5;
gth_sl6 = e1*e2*e3*e4*e5*e6*g0sl6;
gth_sl7 = e1*e2*e3*e4*e5*e6*e7*g0sl7;

%%
% Find elements of body Jacobians via p. 168
fprintf(1, "Find elements of body Jacobians\n");
tw1_dag_sl1 = adinv(e1*g0sl1)*tw1;

tw1_dag_sl2 = adinv(e1*e2*g0sl2)*tw1;
tw2_dag_sl2 = adinv(e2*g0sl2)*tw2;

tw1_dag_sl3 = adinv(e1*e2*e3*g0sl3)*tw1;
tw2_dag_sl3 = adinv(e2*e3*g0sl3)*tw2;
tw3_dag_sl3 = adinv(e3*g0sl3)*tw3;

tw1_dag_sl4 = adinv(e1*e2*e3*e4*g0sl4)*tw1;
tw2_dag_sl4 = adinv(e2*e3*e4*g0sl4)*tw2;
tw3_dag_sl4 = adinv(e3*e4*g0sl4)*tw3;
tw4_dag_sl4 = adinv(e4*g0sl4)*tw4;

tw1_dag_sl5 = adinv(e1*e2*e3*e4*e5*g0sl5)*tw1;
tw2_dag_sl5 = adinv(e2*e3*e4*e5*g0sl5)*tw2;
tw3_dag_sl5 = adinv(e3*e4*e5*g0sl5)*tw3;
tw4_dag_sl5 = adinv(e4*e5*g0sl5)*tw4;
tw5_dag_sl5 = adinv(e5*g0sl5)*tw5;

tw1_dag_sl6 = adinv(e1*e2*e3*e4*e5*e6*g0sl6)*tw1;
tw2_dag_sl6 = adinv(e2*e3*e4*e5*e6*g0sl6)*tw2;
tw3_dag_sl6 = adinv(e3*e4*e5*e6*g0sl6)*tw3;
tw4_dag_sl6 = adinv(e4*e5*e6*g0sl6)*tw4;
tw5_dag_sl6 = adinv(e5*e6*g0sl6)*tw5;
tw6_dag_sl6 = adinv(e6*g0sl6)*tw6;
%%
tw1_dag_sl7 = adinv(e1*e2*e3*e4*e5*e6*e7*g0sl7)*tw1;
%%
tw2_dag_sl7 = adinv(e2*e3*e4*e5*e6*e7*g0sl7)*tw2;
tw3_dag_sl7 = adinv(e3*e4*e5*e6*e7*g0sl7)*tw3;
tw4_dag_sl7 = adinv(e4*e5*e6*e7*g0sl7)*tw4;
tw5_dag_sl7 = adinv(e5*e6*e7*g0sl7)*tw5;
tw6_dag_sl7 = adinv(e6*e7*g0sl7)*tw6;
tw7_dag_sl7 = adinv(e7*g0sl7)*tw7;

%%
% Calculate Jacobians
fprintf(1, "Calculate Jacobians\n");
J1(th1, th2, th3, th4, th5, th6, th7) = simplify([tw1_dag_sl1, zeros(6,6)]);
J2(th1, th2, th3, th4, th5, th6, th7) = simplify([tw1_dag_sl2, tw2_dag_sl2, zeros(6,5)]);
J3(th1, th2, th3, th4, th5, th6, th7) = simplify([tw1_dag_sl3, tw2_dag_sl3,...
    tw3_dag_sl3, zeros(6,4)]);
J4(th1, th2, th3, th4, th5, th6, th7) = simplify([tw1_dag_sl4, tw2_dag_sl4,...
    tw3_dag_sl4, tw4_dag_sl4, zeros(6,3)]);
J5(th1, th2, th3, th4, th5, th6, th7) = simplify([tw1_dag_sl5, tw2_dag_sl5, ...
    tw3_dag_sl5, tw4_dag_sl5, tw5_dag_sl5, zeros(6,2)]);
J6(th1, th2, th3, th4, th5, th6, th7) = simplify([tw1_dag_sl6, tw2_dag_sl6, ...
    tw3_dag_sl6, tw4_dag_sl6, tw5_dag_sl6, tw6_dag_sl6, zeros(6,1)]);
J7(th1, th2, th3, th4, th5, th6, th7) = simplify([tw1_dag_sl7, tw2_dag_sl7, ...
    tw3_dag_sl7, tw4_dag_sl7, tw5_dag_sl7, tw6_dag_sl7, tw7_dag_sl7]);

%Link inertia matrices, via p 172
fprintf(1, "Link inertia matrices\n");
M1 = [m1*eye(3),zeros(3);
    zeros(3), diag(I1)];
M2 = [m2*eye(3),zeros(3);
    zeros(3), diag(I2)];
M3 = [m3*eye(3),zeros(3);
    zeros(3), diag(I3)];
M4 = [m4*eye(3),zeros(3);
    zeros(3), diag(I4)];
M5 = [m5*eye(3),zeros(3);
    zeros(3), diag(I5)];
M6 = [m6*eye(3),zeros(3);
    zeros(3), diag(I6)];
M7 = [m7*eye(3),zeros(3);
    zeros(3), diag(I7)];

% Inertia matrix for the systems, via p 173
fprintf(1, "Inertia matrix for the systems\n");
M = J1'*M1*J1 + J2'*M2*J2 + J3'*M3*J3 + J4'*M4*J4 + J5'*M5*J5 + J6'*M6*J6 + J7'*M7*J7;
fprintf(1, "Mfun\n");
matlabFunction(M,'File','Mfun', 'Optimize', false);
% 

%% Coriolis
fprintf(1, "Coriolis\n");
Mnew = M(th1, th2, th3, th4, th5, th6, th7); % to index as M(i,j)
C = coriolis_barett(Mnew, th1, th2, th3, th4, th5, th6, th7, dth1, dth2,...
    dth3, dth4, dth5, dth6, dth7);
fprintf(1, "Cfun\n");
matlabFunction(C,'File','Cfun', 'Optimize', false);

%%
fprintf(1, "N\n");
N = dV_barrett(th1, th2, th3, th4, th5, th6, th7, m1, m2, m3, m4, m5, m6,...
    m7, gth_sl1, gth_sl2, gth_sl3, gth_sl4, gth_sl5, gth_sl6, gth_sl7);
fprintf(1, "Nfun\n");
matlabFunction(N,'File','Nfun', 'Optimize', false);


% % Calculate inertia matrices
% % iM0 = Q0*diag(I0)*Q0';
% % iM1 = Q1*diag(I1)*Q1';
% % iM2 = Q2*diag(I2)*Q2';
% % iM3 = Q3*diag(I3)*Q3';
% % iM4 = Q4*diag(I4)*Q4';
% % iM5 = Q5*diag(I5)*Q5';
% % iM6 = Q6*diag(I6)*Q6';
% % iM7 = Q7*diag(I7)*Q7';



