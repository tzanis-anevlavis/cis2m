%% Example 1 - "An enhanced hierarchy for (robust) controlled invariant sets"

%% Authors: Tzanis Anevlavis
% Copyright (C) 2019, Tzanis Anevlavis
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
% This script generates an example using the dynamical model of a truck
% with N trailers
%
% This script makes use of the Multi-Parametric Toolbox 3.0:
% M. Herceg, M. Kvasnica, C. Jones, and M. Morari,
% ``Multi-Parametric Toolbox 3.0,'' in Proc. of the European Control
% Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510,
% http://control.ee.ethz.ch/mpt.

close all
clc
clear

%% Model parameters
kd = 4600;  % stiffness
ks = 4500;  % damper coefficient
m0 = 500;   % weight of truck
m = 1000;   % weight of trailer
T = 0.4;    % sampling rate

N = 1;      % number of trailers
disp(N)

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
A = Ao; B = Bo; G = Go; F = Fo;

%% Controlled Invariant Set in Two Moves - Hierarchy
% Start parallel pool if it is not running already:
% if (isempty(gcp('nocreate')))
%     parpool;
% end

% Maximum length of loops:
Loops = [2 3];
Lmax = length(Loops);

% Results matrices:
methods = {'CDC20a','CDC20b'};
numMethods = length(methods);

for mm = 1:numMethods
    Times{mm} = zeros(Lmax,1);
    Volumes{mm} = zeros(Lmax,1);
end

for mm = 1:numMethods  % for each method
    method = methods(mm);
    disp(method)
    
    for l = 1:Lmax      % loops
        L = Loops(l);
        disp(L);

        t = tic;
        cisMat = computeCIS(A,B,G,F,[],[],method,L);
        Times{mm}(l) = toc(t);
        disp(Times{mm}(l));
        cis{mm}(l) = Polyhedron('H', cisMat);
        
        %         % Check if result is numerically invariant:
        %         guard = isInvariant(cis{mm}(l),A,B);
        %         if (~guard)
        %             warning('result not numerically invariant');
        %         end
        
        %         % Check if result is contained by the safe set:
        %         D = Polyhedron('H',[G F]);
        %         guard2 = (cis{mm}(l) <= D);
        %         if (~guard2)
        %             warning('out of safe set');
        %         end
        
        
        % For N>2, volume computation is unstable.
        if (N<3)
            Volumes{mm}(l) = cis{mm}(l).volume;
            disp(Volumes{mm}(l))
        end
    end
    
end