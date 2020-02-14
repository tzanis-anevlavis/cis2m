%% Remark 3 - "A simple hierarchy for computing controlled invariant sets"
% Displays values for Section 5 - Remark 3 - Table 4 - n = 3.

%% Authors: Tzanis Anevlavis, Paulo Tabuada
% Copyright (C) 2019, Tzanis Anevlavis, Paulo Tabuada
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
%
% This code is part of the implementation of the algorithm proposed in:
% Tzanis Anevlavis and Paulo Tabuada, "A simple hierarchy for computing
% controlled invariant sets", Submitted to 23rd ACM International 
% Conference on Hybrid Systems: Computation and Control (HSCC'20),
% and is publicly available at: https://github.com/janis10/cis2m-hierarchy
%
% For any comments contact Tzanis Anevlavis @ janis10@ucla.edu.
%
%
%
%
%% Description:
% This script generates an example using the dynamical model of a truck
% with N trailers, N = 1.
%
% The main code for formulating and solving the optimization problem is
% found on:
% https://homepages.laas.fr/henrion/ (under Software -> ROA)
% to compute an outer approximation of the MCIS for a continuous-time
% linear system.


%% Exammple Setup:
close all
clc
clear

% Model parameters:
kd = 4600;  % stiffness
ks = 4500;  % damper coefficient
m0 = 500;   % weight of truck
m = 1000;   % weight of trailer
T = 0.4;    % sampling rate

%% N = 1 <-> n = 3
N = 1;      % number of trailers
% Continuous-time system - x = [d1 .. dN v0 v1 .. vN] :
tmpA12 = [eye(N) zeros(N,1)]+[zeros(N,1) -eye(N)];
tmpA1 = [zeros(N,N) tmpA12];
tmpA2 = [-ks/m0 zeros(1,N-1) -kd/m0 kd/m0 zeros(1,N-1)];
tmpA31 = ks/m*([eye(N-1) zeros(N-1,1)]+[zeros(N-1,1) -eye(N-1)]);
tmpA32 = kd/m*([eye(N-1) zeros(N-1,1) zeros(N-1,1)]+[zeros(N-1,1) -2*eye(N-1) zeros(N-1,1)]+[zeros(N-1,1) zeros(N-1,1) eye(N-1)]);
tmpA3 = [tmpA31 tmpA32];
tmpA4 = [zeros(1,N-1) ks/m zeros(1,N-1) kd/m -kd/m];
tmpA = [tmpA1; tmpA2; tmpA3; tmpA4];
tmpB = [zeros(N,1); 1; zeros(N,1)];

% Discrete-time system:
% Use forward Euler methnod to obtain discrete-time system, with sampling time T:
A = T*tmpA + eye(2*N+1);
B = T*tmpB;

% State constraints:
d = 0.5;    % spring elongation constraint
vMax = 10;  % max velocity
vMin = 0;   % -min velocity

%% Outer-approximation of the MCIS:

% Degree of v(x) and w(x)
degs = [4 6 8];
for i = 1:length(degs)
    deg = degs(i);
    % Discounting factor
    alpha = 0.9;
    
    % Variables
    n = size(A,1);
    x = sdpvar(n,1);
    u = sdpvar(1,1);
    
    % Dynamics
    f = A*x+B*u;
    
    % Polynomials w(x) and v(x)
    [w, cw] = polynomial(x,deg);
    [v, cv] = polynomial(x,deg);
    
    % This is the only part of this code that we adapted and affects
    % performance in formulating the SDP to be solved. Therefore, we measure
    % the time of solving only the optimization problem later.
    v1 = replace(v,x,f);
    
    % Ball constraints:
    xb1 = d;
    xb2 = vMax;
    xb3 = vMax;
    gx1 = xb1^2 - x(1)^2;
    gx2 = xb2^2 - x(2)^2;
    gx3 = xb3^2 - x(3)^2;
    ub = 5;
    gu = ub^2 - u^2;
    
    % Lebesgue moments on X
    l = getLebesgueMoments(deg,[-d*ones(1,N) -vMin*ones(1,N+1); d*ones(1,N) vMax*ones(1,N+1)],1);
    % Objective
    obj = cw'*l; % minimization w^T l
    
    % Solver options
    SDPsolver = 'mosek'; % or sedumi
    options = getSolverParams(SDPsolver);
    if (~isempty(options))
        mset(options);
    end
    options.verbose = 0;
    
    % SOS multipliers
    % Formulating the constraints also takes some time, not as much as v1 above
    % though.
    [q0,cq0] = polynomial([x;u],deg);
    [q1,cq1] = polynomial([x;u],deg-2);
    [q2,cq2] = polynomial([x;u],deg-2);
    [q3,cq3] = polynomial([x;u],deg-2);
    [r,cr] = polynomial([x;u],deg-2);
    
    [p0,cp0] = polynomial(x,deg);
    [p1,cp1] = polynomial([x;u],deg-2);
    [p2,cp2] = polynomial([x;u],deg-2);
    [p3,cp3] = polynomial([x;u],deg-2);
    
    [s0,cs0] = polynomial(x,deg);
    [s1,cs1] = polynomial([x;u],deg-2);
    [s2,cs2] = polynomial([x;u],deg-2);
    [s3,cs3] = polynomial([x;u],deg-2);
    
    % Constraints
    con = [];
    con = [con; sos( v - alpha*v1 - q0 - gx1*q1 - gx2*q2 - gx3*q3 - gu*r ); sos(q0); sos(q1); sos(q2); sos(q3); sos(r)];
    con = [con; sos( w - v - 1 - p0 - gx1*p1 - gx2*p2 - gx3*p3 ); sos(p0); sos(p1); sos(p2); sos(p3); ];
    con = [con; sos( w - s0 - gx1*s1 - gx2*s2 - gx3*s3 ); sos(s0); sos(s1); sos(s2); sos(s3); ];
    
    % Solve and time *only* for the optimization problem:
    tic
    [~,~,~,res] = solvesos(con,obj,options,[cw;cv;cq0;cq1;cq2;cq3;cr;cp0;cp1;cp2;cp3;cs0;cs1;cs2;cs3]);
    time = toc;
    
    fprintf('Time to solve the SDP for d = %d, is %f seconds. \n',deg,time);
    % fprintf('Residuals %f \n',res);
    % % Coefficients of w(x) and v(x)
    % cw = double(cw);
    % cv = double(cv);
end