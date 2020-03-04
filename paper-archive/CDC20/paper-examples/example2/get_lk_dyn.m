function [A,B,D,U,P,dyn] = get_lk_dyn(param)
% dynamics of lane keeping for fixed longitudinal speed.
% x(t+1) = A x(t) + B u(t) + D d(t)

% parameters
C_af = param.C_af;
C_ar = param.C_ar;
m = param.m;
u = param.u;
b = param.b;
Iz = param.Iz;
a = param.a;
dt = param.dt;

d_f_max = param.steer_max;
d_f_min = param.steer_min;

rd_min = param.rd_min;
rd_max = param.rd_max;

A = eye(4) + [0 1 u 0;
    0 -(C_af+C_ar)/m/u 0 (b*C_ar - a*C_af)/m/u-u;
    0 0 0 1;
    0 (b*C_ar-a*C_af)/Iz/u 0 -(a^2*C_af + b^2*C_ar)/Iz/u]*dt;
B = [0; C_af/m; 0; a*C_af/Iz]*dt;

K = [0;0;0;0];

D = [0;0;-1;0]*dt;

U = Polyhedron('H', [1 d_f_max;
                     -1 -d_f_min]);
                 
XU = Polyhedron('H', [0 0 0 0 1 d_f_max;
                      0 0 0 0 -1 -d_f_min]);
                 
                 
% constraints on the disturbances
P = Polyhedron('H', [1 rd_max;
                     -1 -rd_min]);
                 
dyn = Dyn(A, K, B, XU, [],[],[],{zeros(4,4)}, {D}, P);    

end