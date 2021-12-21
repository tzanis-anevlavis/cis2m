function [mcisA,mcisb, A_hd, K, P] = implicitclosedformRCIS(Ac,Bc,G_k,F_k,Lambda,Tau,nmax)
%% Authors: T.Anevlavis, Z.Liu, N.Ozay, and P.Tabuada
% Copyright (C) 2021, T.Anevlavis, Z.Liu, N.Ozay, and P.Tabuada
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
% For any comments contact Tzanis Anevlavis @ t.anevlavis@ucla.edu.
%
%
%
%

n = size(Ac,2); m = size(Bc,2);

%% Construct the high-dimensional dynamical system:
K = kron(speye(m),[1 sparse(1,Tau+Lambda-1)]);
P = kron(speye(m),[sparse(Tau+Lambda-1,1) speye(Tau+Lambda-1); sparse(1,Tau) 1 sparse(1,Lambda-1)]);

A_hd = [Ac Bc*K; sparse((Tau+Lambda)*m,n) P];

%% Construct high-dimensional invariant set:

% blkdiag([Gc 0], [G_k{1} 0], .. , [G_k{nmax} 0], .. , [G_k{nmax} 0])
blk_G_k = blkdiag(G_k{:});

A_mats = cell(nmax+(Tau+Lambda),1);
A_mats{1} = speye(size(A_hd));
for t = 2:nmax+(Tau+Lambda)
    A_mats{t} = A_mats{t-1} * A_hd;
end
A_mats = cat(1, A_mats{:});

mcisA = blk_G_k * A_mats;

mcisb = cat(1, F_k{:}); % [F; F_k{1}; .. F_k{nmax}; .. F_k{nmax}]

%% Check if positively invariant:
mcis = Polyhedron('H', [mcisA mcisb]);
if (~isPositivelyInvariant(mcis, A_hd))
    warning("Implicit CIS not positively invariant in Brunovsky space!");
end
