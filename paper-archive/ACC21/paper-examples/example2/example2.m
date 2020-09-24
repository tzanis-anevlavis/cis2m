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
% Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510,
% http://control.ee.ethz.ch/mpt.
%
% This script makes use of the PCIS toolbox:
% https://github.com/pettni/pcis

close all
clear
clc

%% Model parameters:
% Original system:
[param,Safe] = constant_lk4();
[A,B,D,U,P] = get_lk_dyn(param);

% State constraints:
G = Safe.A;
F = Safe.b;

% Input constraints:
Gu = U.A;
Fu = U.b;

%% Compute MCIS using MPT3
tic
system = LTISystem('A',A,'B',B);
mcis = system.invariantSet('X',Safe,'U',U,'maxIterations',300);
timeMCIS = toc;
volumeMCIS = mcis.volume;

%% Controlled Invariant Set in Two Moves - Hierarchy

% Maximum length of loops:
Loops = (2:6);
Lmax = length(Loops);

% Results matrices:
methods = {'ACC21a','ACC21b'};
numMethods = length(methods);

for mm = 1:numMethods
    Times{mm} = zeros(Lmax,1);
    Volumes{mm} = zeros(Lmax,1);
    VolumePercentage{mm} = zeros(Lmax,1);
end

% Start parallel pool if it is not running already:
if (isempty(gcp('nocreate')))
    parpool;
end

for mm = 1:numMethods
    method = methods(mm);
    disp(method)
    
    for l = 1:Lmax
        L = Loops(l);
        disp(L);
        
        t2  = tic;
        cisMat = computeCIS(A,B,G,F,Gu,Fu,method,L);
        Times{mm}(l) = toc(t2);
        disp(Times{mm}(l));
        cis{mm}(l) = Polyhedron('H', cisMat);
        
        guard = isInvariant(cis{mm}(l),A,B);
        if (~guard)
            warning('result not numerically invariant');
        end
        
        D = Polyhedron('H',[G F]);
        guard2 = (cis{mm}(l) <= D);
        if (~guard2)
            warning('out of safe set');
        end
        
        Volumes{mm}(l) = cis{mm}(l).volume;
        disp(Volumes{mm}(l));
    end
    
    VolumePercentage{mm} = Volumes{mm}./volumeMCIS*100;
end
