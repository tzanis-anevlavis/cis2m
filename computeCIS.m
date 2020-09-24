function [cisMat] = computeCIS(A,B,G,F,Gu,Fu,method,L,E,Gw,Fw,verbose)
%% Authors: Tzanis Anevlavis
% Copyright (C) 2020, Tzanis Anevlavis
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
% This is the wrapper of the CIS2M repository.
%
% Takes as input a discete-time linear system and a polyhedral safe set.
% Computes a  (robust) Controlled Invariant Subset of the safe set.
%
% Inputs:   A, B : matrices defining the discrete-time LTI: x+ = Ax + Bu
%           G, F : matrices defining the safe set: {x \in \R^n | G x <= F}
%           Gu,Fu: matrices defining the input constr: {u \in \R | Gu u <= Fu}
%                  if there are no input constraints use Gu = [], Fu = [].
%              E : disturbance matrix if the system is: x+ = Ax + Bu + Ew
%           Gw,Fw: matrices defining the disturb. set: {w \in \R^k | Gw w <= Fw}
%                  if there is no disturbance use E = [], Gu = [], Fu = [].
%
%           method: available methods CDC19, HSCC20, ACC21a, ACC21b.
%                  in the case of disturbance choose ACC21b.
%
%              L : lentgh of the loop for HSCC20, ACC21a, ACC21b (L>0 int.)
%
%           verbose = 0 - no messages; 1 - displays messages.
%
% Output:   cisMat: matrix defining the computed Controlled Invariant Set
%                   cisMat = [cisA cib]
%                   s.t. CIS = {x \in \R^n | cisA x <= cisb}

%% Add support folder to path
% addpath('./support_functions/');

%% Input arguments check
n = size(A,1);	% dimension of the system
inputArgCheck

%% Convert system in Brunovsky normal form space:
[Ac,Bc,Ec,Gc,Pmat] = convert2Bru(A,B,E,G);

%% Input constraints:
% If there are input constraints, we extend the system by one dimension to
% incorporate them.
guard = false;
if (~isempty(Gu))
    guard = true;
    [A,B,E,G,F] = inputExt(A,Pmat,Ac,Bc,Ec,Gc,F,Gu,Fu);
    n = size(A,1);
else
    A = Ac; B = Bc; E = Ec; G = Gc;
end

%% Construct D_k
% This is for ACC21b only.
if (disturbance)
    % Minkowski difference for the closed-form expression.
    G_k = cell(n,1);
    F_k = cell(n,1);
    W = Polyhedron('A',Gw,'b',Fw);
    D = Polyhedron('A',G,'b',F);
    for i = 1:n
        D = D - A^(i-1)*E*W;
        D.minHRep();
        G_k{i} = D.A;
        F_k{i} = D.b;
    end
else
    G_k = cell(n,1);
    F_k = cell(n,1);
    D = Polyhedron('A',G,'b',F);
    D.minHRep();
    for i = 1:n
        G_k{i} = D.A;
        F_k{i} = D.b;
    end
end

%% Controlled Invariant Set in Two Moves
if (strcmp(method,'CDC19'))
    [cisLiftedA,cisLiftedb] = cdc19(A,B,G,F,verbose);
elseif (strcmp(method,'HSCC20'))
    if (L==1)
        [cisLiftedA,cisLiftedb] = cdc19(A,B,G,F,verbose);
    else
        [cisLiftedA,cisLiftedb] = hscc20(A,B,G,F,L,verbose);
    end
elseif (strcmp(method,'ACC21a'))
    [cisLiftedA,cisLiftedb] = acc21a(A,G,F,L,verbose);
elseif (strcmp(method,'ACC21b'))
    [cisLiftedA,cisLiftedb] = acc21b(A,G,F,G_k,F_k,L,verbose);
else
    error('Invalid method.');
end

proj = 0;
% Use custom projection (faster with parallel pool)
if (proj == 1)
    [cisA,cisb] = jProject(cisLiftedA,cisLiftedb,n,verbose);
    % If we had input constraints, we project from the extended space to the
    % original space, i.e., eliminate input u.
    if (guard)
        [cisA, cisb] = jProject(cisA,cisb,n-1,verbose);
    end
    % Use MPT3 and let it decide.
else
    P = Polyhedron('A',cisLiftedA,'b',cisLiftedb);
    P.minHRep;
%     Px = P.projection(1:n);
    Px = P.projection(1:n,'ifourier');  % seems to be better than mplp for
                                        % many cases.

    % If we had input constraints, we project from the extended space to the
    % original space, i.e., eliminate input u.
    if (guard)
        Px = Px.projection(1:(n-1));
    end
    
    cisA = Px.A;
    cisb = Px.b;
end

cisMat = [cisA cisb];

% Return to original coordinates:
cisMat = [cisMat(:,1:end-1)*Pmat cisMat(:,end)];


%% Remove support folder from path
% rmpath('./support_functions/');