function [mcisA,mcisb] = mcisCF(Ac,Bc,Gc,F,verbose)

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
% Tzanis Anevlavis and Paulo Tabuada, "Computing controlled invariant sets
% in two moves", in 2019 IEEE Conference on Decision and Control,
% and is publicly available at: https://github.com/janis10/cis2m
%
% For any comments contact Tzanis Anevlavis @ janis10@ucla.edu.
%
%
%
%
%% Description:
% This function takes as input a system in Brunovsky normal form 
% and a polyhedron in the corresponding coordinates, and returns the
% Maximal Controlled Invariant Set (MCIS) (in the higher dimensional space)
%
% System:                   x+ = Ac x + Bc u
% Polyhedral constraint:    D = {x \in \R^n | Gc x <= F}
% 
% Returns: matrices mcisA, mcisb such that
%                   MCIS = {(x,f)| mcisA (x,f)^T <= mcisb}
%
%           verbose = 0 - no messages; 1 - displays messages.

n = size(Ac,1); % number of original variables
m = size(Gc,1); % number of original inequalities

%% Construct MCIS in High Dimension
if (verbose)
    disp('Lifting problem to compute controlled invariant set in closed-form . . .')
end

% Construct Ghat0 and Hhat0:
tmpG = zeros(m*n,n);
Hhat0 = zeros(m,n+m*n);
for j = 1:m
    % For Ghat0:
    tmpG((j-1)*n+1:j*n,:) = diag(Gc(j,:));
    
    % Hhat0 here --faster than doing sparse.
    for i = 1:n
        Hhat0(j,j*n+i) = 1;
    end
end
tmpG = sparse(tmpG);
Ghat0 = [tmpG -speye(m*n)];
Hhat0 = sparse(Hhat0);

% "Index" matrices:
Gsign = Gc';
Gsign = Gsign(:);
P = zeros(n*m,1);
T = zeros(n*m,1);
for i = 1:n*m
    if (Gsign(i)>0)
        P(i) = -1/Gsign(i);
    elseif (Gsign(i)<0)
        T(i) = 1/Gsign(i);
    end
end
Pidx = find(P~=0);
Tidx = find(T~=0);

% Construct the combinatorial constraints matrix
sP = size(Pidx,1);
sT = size(Tidx,1);
hGamma = sP*sT;

gammaCell = cell(sP,1);
% pre-allocate
for p = 1:sP
    gammaCell{p} = sparse(n+m*n,sT);
end
parfor p = 1:sP
    for t = 1:sT
        gammaCell{p}(n+Pidx(p),t) = P(Pidx(p));
        gammaCell{p}(n+Tidx(t),t) = T(Tidx(t));
    end
end
Gamma = [gammaCell{:}]';
    
% Costruct the passing-through-dynamics-constraints matrix Gtilde
AcSP = sparse(Ac);
Ahat = [AcSP sparse(n,n*m); sparse(n*m,n) speye(n*m,n*m)];
BcSP = sparse(Bc);
Bhat = [BcSP; sparse(n*m,1)];
% Ridx = indices that involve 'u'.
% tempA1 = Ghat0 Ahat^{1} till Ghat0 Ahat^{n} at each step 'i' without
% indices in Ridx1.
% tempB1 = Ghat0 Ac^{0} Bc till Ghat0 Ac^{n-1} Bc at each step 'i' without
% indices in Ridx1.
% Gtilde = concatinate vertically the above resulting inequalities.

gtildeCell = cell(n,1);
tempA1 = Ghat0*Ahat;
tempB1 = Ghat0*Bhat;
% gotta fix this: it keeps by mistake rows of the form: -\lambda_ji < 0
% can do (n-1) cause at n it's all zero.
% fixed and updated Git!
for i = 1:n
    Ridx = (tempB1~=0);
    gtildeCell{i} = tempA1;
    gtildeCell{i}(Ridx,:) = [];
    gtildeCell{i} = gtildeCell{i}';
    
    tempA1(Ridx,:)=[];
    tempB1 = tempA1*Bhat;
    tempA1 = tempA1*Ahat;
end
Gtilde = [gtildeCell{:}]';
Ftilde = sparse(size(Gtilde,1),1);

% MCIS in closed-form:
% [Sum of f vars; Original lift; Constraints by passing through dynamics;
% Constraints on f vars from FME]
mcisA = [Hhat0; Ghat0; Gtilde; Gamma];

F = sparse(F);

mcisb = [F; sparse(m*n,1); Ftilde; sparse(hGamma,1)];

if (verbose)
    disp('..computed!')
end