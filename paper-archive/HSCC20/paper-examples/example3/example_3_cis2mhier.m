%% Example 3 - "A simple hierarchy for computing controlled invariant sets"
% Values for Section 5 - Example 3 - Table 2 - Algorithm 1:
% Times
% Volumes
% Values for Section 5 - Example 3 - Table 2 - Algorithm of [23]:
% VolElls

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
% clc
clear

% Loop length:
L = 2;

Nmax = 4;   % max number of trailers
% Results matrices:
Times = zeros(Nmax,1);         % stores running times
Volumes = zeros(Nmax,1);       % stores computed volumes

%%%%%%%%%%%%%%%%%%%%%%%% Model parameters %%%%%%%%%%%%%%%%%%%%%%%%
kd = 4600;  % stiffness
ks = 4500;  % damper coefficient
m0 = 500;   % weight of truck
m = 1000;   % weight of trailer
T = 0.4;    % sampling rate

for N = 1:4      % number of trailers
    disp(N);
    
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
    
    %%%%%%%%%%% Controlled Invariant Set in Two Moves - Hierarchy %%%%%%%%%%%
    % Start parallel pool if it is not running already:
    if (isempty(gcp('nocreate')))
        parpool;
    end
    
    method = 'HSCC20';
    
    tic
    cisMat = computeCIS(A,B,G,F,[],[],method,L);
    Times(N) = toc;
	cis(N) = Polyhedron('H', cisMat);
    
    if (N == 3) % we noticed that MPT3's volume function was returning an
        % error when trying to compute the volume in this case, but
        % works fine if we decrease precision of decimal digits of
        % the vertex representation.
        cisVert = cis(N).V;
        cisVert = cisVert.*1e4;
        cisVert = round(cisVert);
        cisVert = cisVert./1e4;
        cis(N) = Polyhedron('V',cisVert);
    end
    
    if (N<4)    % For N=4, volume computation exceeded 5 hours.
        Volumes(N) = cis(N).volume;
        disp(Volumes(N))
    end
end

%% Compute the volumes of controlled invariant ellipsoids from
% B.Legat et al "Computingcontrolled invariant sets for hybrid systems with
% applications to model-predictive control." IFAC ADHS 2018.

% An ellipsoid is defined as {x\in\R^n | x^T P x <=1}.
% The values for P matrices below were computed in the notebook
% "example_3_ellips_times.ipynb", which implements the algorithm of the
% above paper, and is provided in folder "./paper-examples/ellipsoidalCIS".

% For N = 1 <-> n = 3:
P1 = [0.249984 -0.224672 0.0790755; -0.224672 0.992715 0.77826; 0.0790755 0.77826 1.22903];
singVal1 = sqrt(eig(P1));
volEll1 = pi^(3/2)*prod(singVal1)/gamma(3/2+1);
% For N = 2 <-> n = 5:
P2 = [0.24988 -0.0636891 -0.0964877 0.163941 0.0972664; -0.0636891 0.249897 -0.0945175 -0.183021 0.0832756; -0.0964877 -0.0945175 0.686088 0.452582 0.183919; 0.163941 -0.183021 0.452582 0.741792 0.577272; 0.0972664 0.0832756 0.183919 0.577272 0.881588];
singVal2 = sqrt(eig(P2));
volEll2 = pi^(5/2)*prod(singVal2)/gamma(5/2+1);
% For N = 3 <-> n = 7:
P3 = [0.247311 -0.0197818 -0.00201309 -0.111419 0.120498 0.0510998 0.00487735; -0.0197818 0.247695 -0.0804444 -0.0556228 -0.0962673 0.163089 0.0840975; -0.00201309 -0.0804444 0.248017 -0.0403077 -0.0981401 -0.195127 0.0776465; -0.111419 -0.0556228 -0.0403077 0.535786 0.389067 0.200226 0.0508651; 0.120498 -0.0962673 -0.0981401 0.389067 0.62922 0.420648 0.165395; 0.0510998 0.163089 -0.195127 0.200226 0.420648 0.665933 0.468254; 0.00487735 0.0840975 0.0776465 0.0508651 0.165395 0.468254 0.727407];
singVal3 = sqrt(eig(P3));
volEll3 = pi^(7/2)*prod(singVal3)/gamma(7/2+1);

VolElls = [volEll1 volEll2 volEll3];