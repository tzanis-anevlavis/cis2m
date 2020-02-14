%% Example 3 - "A simple hierarchy for computing controlled invariant sets"
% Values for Section 5 - Example 3 - Table 2 - MPT:
% TimeMCIS
% VolumeMCIS

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
% with N trailers
%
% This script makes use of the Multi-Parametric Toolbox 3.0:
% M. Herceg, M. Kvasnica, C. Jones, and M. Morari,
% ``Multi-Parametric Toolbox 3.0,'' in Proc. of the European Control
% Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510,
% http://control.ee.ethz.ch/ mpt.

%% Exammple Setup:
close all
clc
clear

Nmax = 4;   % max number of trailers
% Results matrices:
TimeMCIS = zeros(Nmax,1);           % stores running times
VolumeMCIS = zeros(Nmax,1);         % stores volumes
%%%%%%%%%%%%%%%%%%%%%%%% Model parameters %%%%%%%%%%%%%%%%%%%%%%%%
kd = 4600;  % stiffness
ks = 4500;  % damper coefficient
m0 = 500;   % weight of truck
m = 1000;   % weight of trailer
T = 0.4;    % sampling rate

for N = 1:Nmax      % number of trailers
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
    Ao = T*tmpA + eye(2*N+1);
    Bo = T*tmpB;
    
    % State constraints:
    d = 0.5;    % spring elongation constraint
    vMax = 10;  % max velocity
    vMin = 0;   % -min velocity
    
    Go =[eye(N) zeros(N,N+1); -eye(N) zeros(N,N+1); zeros(N+1,N) eye(N+1); zeros(N+1,N) -eye(N+1)];
    Fo = [d*ones(N,1); d*ones(N,1); vMax*ones(N+1,1); vMin*ones(N+1,1)];
    D = Polyhedron('H',[Go Fo]);
    
    % Final matrices:
    A = Ao; B = Bo;
    
    %%%%%%%%%%%%%%%%%%%%%%%% Compute MCIS using MPT3 %%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    system = LTISystem('A',A,'B',B);
    mcis = system.invariantSet('X',D,'maxIterations',50);
    TimeMCIS(N) = toc;
    if (N<3)    % For N = 3,4 the algorithm does not converge and we do not
                % compute the volumes.
        VolumeMCIS(N) = mcis.volume;
        disp(VolumeMCIS(N))
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end