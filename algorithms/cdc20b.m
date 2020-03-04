function [mcisA,mcisb] = cdc20b(Ac,Gc,F,G_k,F_k,L,verbose)
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
%% Description:
% Computes a lifted controlled invariant set using Algorithm 2, proposed in
% T.Anevlavis, Z.Liu, N.Ozay, and P.Tabuada, "An enhanced hierarchy for
% (robust) controlled invariant sets", Submitted to 2020 IEEE Conference on
% Decision and Control (CDC'20).


if (verbose)
    disp('Lifting problem to compute controlled invariant set in closed-form . . .')
end

n = size(Ac,1);
k = size(Gc,1);

dbar = floor(n/L);
S = [zeros(L-1,1) eye(L-1); 1 zeros(1,L-1)];
I = repmat(eye(L),dbar,1);
P = [I; eye(n-dbar*L) zeros(n-dbar*L,(1+dbar)*L-n)];

%% Step 0 constraints:
% Gc x \leq F
G0 = [Gc zeros(k,L)];
F0 = F;

%% Step 1 to n-1 constraints:
G1n1 = [];
F1n1 = [];

for s = 1:n-1
    sbar = min(s,n);
    T = [[zeros(n-sbar,sbar); eye(sbar)] zeros(n,n-sbar) ];
    
    G1n1 = [G1n1; G_k{s}*Ac^s G_k{s}*T*P];
    F1n1 = [F1n1; F_k{s}];
end

%% Step n to n+L-1 constraints:
GnL1 = [];
FnL1 = [];

for s = n:n+L-1
    sbar = min(s,n);
    T = [[zeros(n-sbar,sbar); eye(sbar)] zeros(n,n-sbar) ];
    
    GnL1 = [GnL1; zeros(size(G_k{n},1),n) G_k{n}*T*P*S^(s-1)];
    FnL1 = [FnL1; F_k{n}];
end

%% Construct the final set:
finalSet = [G0 F0; G1n1 F1n1; GnL1 FnL1];

mcisA = finalSet(:,1:end-1);
mcisb = finalSet(:,end);

if (verbose)
    disp('..computed!')
end