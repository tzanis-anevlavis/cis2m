%% Example 2 - "An enhanced hierarchy for (robust) controlled invariant sets"

%% Authors: Tzanis Anevlavis, Zexiang Liu
% Copyright (C) 2020, Tzanis Anevlavis, Zexiang Liu
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
%
% For any comments contact Tzanis Anevlavis @ janis10@ucla.edu.
%
%
%
%
%% Description:
% This script generates a 4D lane keeping example.
%
% This script makes use of the Multi-Parametric Toolbox 3.0:
% M. Herceg, M. Kvasnica, C. Jones, and M. Morari,
% ``Multi-Parametric Toolbox 3.0,'' in Proc. of the European Control
% Conference, Z�rich, Switzerland, July 17-19 2013, pp. 502-510,
% http://control.ee.ethz.ch/mpt.
%
% This script makes use of the PCIS toolbox:
% https://github.com/pettni/pcis

close all
clear
clc
% Relative to this file:
addpath('../../../legacy_code/');

%% Generate system dynamics
[param,Safe] = constant_lk4();  % Choose disturbance range within. 
[A,B,E,U,P,dyn] = get_lk_dyn(param);

Gx = Safe.A;
Fx = Safe.b;
Gu = U.A;
Fu = U.b;
Gw = P.A;
Fw = P.b;

%% Compute MCIS using MPT3
t = tic;
mcis = dyn.win_always(Safe,0,0,1);
timeMCIS = toc(t);
disp(timeMCIS)
volumeMCIS = mcis.volume;
disp(volumeMCIS)

%% Controlled Invariant Set in Two Moves - Hierarchy

% Start parallel pool if it is not running already:
if (isempty(gcp('nocreate')))
    parpool;
end

% Length of loops:
% Loops = (1:6);
Loops = 4;
Lmax = length(Loops);

% Results matrices:
Times = zeros(Lmax,1);
Volumes = zeros(Lmax,1);
VolumePercentage = zeros(Lmax,1);

method = 'ACC21b';
disp(method)
for l=1:Lmax
    L = Loops(l);
    disp(L)
    
    t = tic;
    cisMat = computeCIS(A,B,Gx,Fx,Gu,Fu,E,Gw,Fw,method,L);
    Times(l) = toc(t);
    
    cis(l) = Polyhedron('H',cisMat);
    Volumes(l) = cis(l).volume;
end

disp(Times)

VolumePercentage = Volumes{mm}./volumeMCIS*100;

disp(VolumePercentage)