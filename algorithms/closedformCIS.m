function [mcisA,mcisb] = closedformCIS(Ac,Bc,Gc,F,G_k,F_k,L,nmax,verbose)
%% Authors: T.Anevlavis, Z.Liu, N.Ozay, and P.Tabuada (alphabetically)
% Copyright (C) 2020, T.Anevlavis, Z.Liu, N.Ozay, and P.Tabuada
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
% This code is part of the Controlled Invariance in 2 Moves repository 
% (CIS2M), and is publicly available at: https://github.com/janis10/cis2m .
%
% For any comments contact Tzanis Anevlavis @ janis10@ucla.edu.
%
%
%
%


if (verbose)
    disp('Lifting problem to compute controlled invariant set in closed-form . . .')
end

n = size(Ac,2);
m = size(Bc,2);

%% Construct the high-dimensional dynamical system:
Ti = [1 zeros(1,L-1)];
T = [];
Pi = [zeros(L-1,1) eye(L-1); 1 zeros(1,L-1)];
P = [];
for i = 1:m
    T = blkdiag(T,Ti);
    P = blkdiag(P,Pi);
end

% [ z+ ]        [   Abru    Bbru    0   ] [z]
% [ u+ ]  =     [   0       0       T   ] [u]
% [ v+ ]        [   0       0       P   ] [v]
A_hd = [Ac             Bc*T;
        zeros(L*m,n)    P   ];

%% Construct high-dimensional invariant set:
% Initial:
mcisA = [Gc zeros(size(Gc,1),m*L)];
mcisb = F;
A_curr = A_hd;
% Constrained reachability:
for t = 1:nmax+L-1
    tbar = min(t,nmax);
    mcisA = [mcisA; [G_k{tbar} zeros(size(G_k{tbar},1),m*L)]*A_curr];
    mcisb = [mcisb; F_k{tbar}];
    A_curr = A_hd * A_curr;
end

if (verbose)
    disp('..computed!')
end
