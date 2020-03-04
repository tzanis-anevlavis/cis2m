%% Example 1 & 2 - "A simple hierarchy for computing controlled invariant sets"
% Plots for Section 2 - Example 1.
% Values of Algorithm 1 for Section 5 - Example 2 - Table 1:
% Times
% VolumePercentage

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
% and is publicly available at: https://github.com/janis10/cis2m
%
% For any comments contact Tzanis Anevlavis @ janis10@ucla.edu.
%
%
%
%
%% Description:
% This script generates an example in \R^2, with constrained input u.
% Constraints are incorporated by extending the original system by one
% dimension obtaining the state y = (x,u), with unconstrained input v.
% Therefore, the resulting system for our algorithm is in \R^3.
%
% Original system: x+ = Ao x + Bo
% State constraints: Go x <= Fo
% Input constraints: u \in [-umin,umax].
%
% This script makes use of the Multi-Parametric Toolbox 3.0:
% M. Herceg, M. Kvasnica, C. Jones, and M. Morari,
% ``Multi-Parametric Toolbox 3.0,'' in Proc. of the European Control
% Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510,
% http://control.ee.ethz.ch/ mpt.

%% Exammple Setup:
close all
clear
clc

% Original system:
A = [1.5 1; 0 1]; B = [0.5; 0.25];
% State constraints:
G = [   0.9147   -0.5402;
    0.2005    0.6213;
    -0.8193    0.9769;
    -0.4895   -0.8200;
    0.7171   -0.3581;
    0.8221    0.0228;
    0.3993   -0.8788];
F = [   0.5566;
    0.8300;
    0.7890;
    0.3178;
    0.4522;
    0.7522;
    0.1099];
D = Polyhedron('H',[G F]);
% Input constraints:
umin = -1;
umax = 1;
Gu = [1; -1]; Fu = [umax; -umin];
U = Polyhedron('H',[Gu Fu]);

%%%%%%%%%%%%%%%%%%%%%%%% Compute MCIS using MPT3 %%%%%%%%%%%%%%%%%%%%%%%%%
tic
system = LTISystem('A',A,'B',B);
mcis = system.invariantSet('X',D,'U',U,'maxIterations',300);
timeMCIS = toc;
volumeMCIS = mcis.volume;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Controlled Invariant Set in Two Moves - Hierarchy %%%%%%%%%%%

% Maximum length of loops:
Lmax = 6;
Times = zeros(Lmax,1);
Volumes = zeros(Lmax,1);

% Start parallel pool if it is not running already:
if (isempty(gcp('nocreate')))
    parpool;
end

for L = 1:Lmax
    disp(L);
    
    method = 'HSCC20';
    tic
    cisMat = computeCIS(A,B,G,F,Gu,Fu,method,L);
    Times(L) = toc;
    cis(L) = Polyhedron('H', cisMat);
    
    Volumes(L) = cis(L).volume;
end
VolumePercentage = Volumes./volumeMCIS*100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Colors:
greenShades = zeros(4,3);
greenShades(1,:) = 1/255*[0,255,128];
greenShades(2,:) = 1/255*[0,204,102];
greenShades(3,:) = 1/255*[0,153,76];
greenShades(4,:) = 1/255*[0,102,51];

figure; plot(D,'color','blue',mcis, 'color', 'lightgray', cis(2), 'color', greenShades(1,:), cis(1), 'color', 'white')
axis([-0.9 1 -0.3 1.2])
figure; plot(mcis, 'color', 'lightgray', cis(3), 'color', greenShades(2,:), cis(2), 'color', greenShades(1,:), cis(1), 'color', 'white')
axis([-0.9 1 -0.3 1.2])
figure; plot(mcis, 'color', 'lightgray', cis(4), 'color', greenShades(3,:), cis(3), 'color', greenShades(2,:), cis(2), 'color', greenShades(1,:), cis(1), 'color', 'white')
axis([-0.9 1 -0.3 1.2])
figure; plot(mcis, 'color', 'lightgray', cis(5), 'color', greenShades(4,:), cis(4), 'color', greenShades(3,:), cis(3), 'color', greenShades(2,:), cis(2), 'color', greenShades(1,:), cis(1), 'color', 'white')
axis([-0.9 1 -0.3 1.2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save and Export %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G0 = mcis.A;    F0 = mcis.b;
% G1 = cis(1).A;  F1 = cis(1).b;
% G2 = cis(2).A;  F2 = cis(2).b;
% G3 = cis(3).A;  F3 = cis(3).b;
% G4 = cis(4).A;  F4 = cis(4).b;
% G5 = cis(5).A;  F5 = cis(5).b;
% filename = 'mat2julia.mat';
% save(filename,'G0','F0','G1','F1','G2','F2','G3','F3','G4','F4','G5','F5')