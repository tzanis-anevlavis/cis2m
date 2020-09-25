function [param,Safe] = constant_lk4()
% safety constraints      
y_max = 0.9;
y_min = -0.9;
v_max = 1.2;
v_min = -1.2;
Phi_max = 0.05;
Phi_min = -0.05;
dPhi_max = 0.3;
dPhi_min = -0.3;

Safe = Polyhedron("H",...
    [1 0 0 0 y_max
    -1 0 0 0 -y_min
    0 1 0 0 v_max
    0 -1 0 0 -v_min
    0 0 1 0 Phi_max
    0 0 -1 0 -Phi_min
    0 0 0 1 dPhi_max
    0 0 0 -1 -dPhi_min]);

% vehicle parameters
param.C_af = 133000;
param.C_ar = 98800;
param.m = 1650;
param.u = 30;
param.b = 1.59;
param.Iz = 2315;
param.a = 1.11;
param.steer_max = pi/2; 
param.steer_min = -pi/2;

param.dt = 0.1;

% disturbance range
param.rd_min = -0.02;
param.rd_max = 0.02;