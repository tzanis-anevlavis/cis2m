function [cisA, cisb] = computeCIS(A,B,Gx,Fx,Gu,Fu,E,Gw,Fw,method,L,verbose)
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
% Inputs:   A,B,E : matrices defining the discrete-time LTI: x+ = Ax + Bu + Ew
%           Gx,Fx: define the safe set: {x \in \R^n | Gx x <= Fx}
%           Gu,Fu: define the input constr: {u \in \R | Gu u <= Fu}
%                  if there are no input constraints use Gu = [], Fu = [].
%           Gw,Fw: define the disturbance set: {w \in \R^k | Gw w <= Fw}
%                  if there is no disturbance use E = [], Gw = [], Fw = [].
%
%           method: legacy options for reference to our papers:
%                       CDC19, HSCC20, ACC21a, ACC21b.
%
%              L : denotes the level of the hierarchy (L-invariance).
%
%           verbose: 0 - no messages; 1 - displays messages.
%
% Output:   [cisA, cib]: matrix defining the computed Controlled Invariant Set
%                   s.t. CIS = {x \in \R^n | cisA x <= cisb}

%% Input arguments check
% If no verbose, choose silent mode.
if (~exist('verbose','var'))
    verbose = 0;
end
% Check state constraints and system:
if (~exist('Gx','var') && ~exist('Fx','var'))
    error('State constraints not specified.')
elseif (exist('Gx','var') && ~exist('Fx','var'))
    error('State constraints incomplete: Matrix Gx given, but not vector Fx.')
elseif (~exist('Gx','var') && exist('Fx','var'))
    error('State constraints incomplete: Vector Fx given, but not matrix Gx.')
end
if (size(Gx,1)~=size(Fx,1))
    error('Rows of Gx and Fx do not match.')
elseif (size(Gx,2)~=size(A,2))
    error('Columns of A (number of states) and Gx do not match.')
elseif (size(A,1)~=size(B,1))
    error('Rows of A and B do not match.')
end
% Check input constraints:
if (~exist('Gu','var') && ~exist('Fu','var'))
    Gu = [];    Fu = [];
elseif (exist('Gu','var') && ~exist('Fu','var'))
    error('Input constraints incomplete: Matrix Gu given, but not vector Fu.')
elseif (~exist('Gu','var') && exist('Fu','var'))
    error('Input constraints incomplete: Vector Fu given, but not matrix Gu.')
elseif (~isempty(Gu) && ~isempty(Fu))
    if (size(Gu,1)~=size(Fu,1))
        error('Rows of Gu and Fu do not match.')
    elseif (size(B,2)~=size(Gu,2))
        error('Columns of B (number of inputs) and Gu do not match.')
    end
end
% Disturbance matrix and set. If no disturbance set E=[], Gw=[], Fw=[].
if (~exist('E','var') || isempty(E))
    E = []; Gw = []; Fw = [];
else
    if (~exist('Gw','var') && ~exist('Fw','var'))
        error('Disturbance matrix E is specified, but not disturbance set.')
    elseif (exist('Gw','var') && ~exist('Fw','var'))
        error('Disturbance set incomplete: Matrix Gw given, but not vector Fw.')
    elseif (~exist('Gw','var') && exist('Fu','var'))
        error('Disturbance set incomplete: Vector Fw given, but not matrix Gw.')
    elseif (~isempty(Gw) && ~isempty(Fw))
        if (size(Gw,1)~=size(Fw,1))
            error('Rows of Gw and Fw do not match.')
        elseif (size(E,2)~=size(Gw,2))
            error('Columns of E and Gw do not match.')
        else
            W = Polyhedron('H',[Gw Fw]);
            if (~W.isBounded)
                error('Disturbance set is unbounded.')
            end
        end
    end
end
% If no method selected, error:
if (~exist('method','var'))
    error('No method selected.'); 
end
% Method compatibility with disturbance case:
if(disturbance)
    if ( (strcmp(method,'CDC19')) || (strcmp(method,'HSCC20')) || (strcmp(method,'ACC21a')) )
        error('The selected method has not been implemented for the case of disturbances. Please select ACC21b.')
    end
end
% Loop selection:
if (~(exist('L','var') && (L>0) && (floor(L)==L) ) )
    disp('Length of the loop (L) must be a positive integer.')
    prompt = 'Choose a value for length of loop (L>0):';
    L = input(prompt);
    while ((L<0) || (L==0) || (floor(L)~=L))
        prompt = 'Input L must be a positive integer:';
        L = input(prompt);
    end
end
if (strcmp(method,'CDC19'))
    L = 1;
end

%% Convert system in Brunovsky normal form space, and extend state space:
[Ac,Bc,Ec,Gc,Fc,Pmat,nmax,~,~,isExtended] = convert2Bru(A,B,E,Gx,Fx,Gu,Fu)
n = size(Ac,2);
m = size(Bc,2);

%% Construct D_k
if (disturbance)
    % Minkowski difference for the closed-form expression.
    G_k = cell(nmax,1);
    F_k = cell(nmax,1);
    W = Polyhedron('A',Gw,'b',Fw);
    D = Polyhedron('A',Gc,'b',Fc);
    for i = 1:nmax
        D = D - Ac^(i-1)*Ec*W;
        D.minHRep();
        G_k{i} = D.A;
        F_k{i} = D.b;
    end
else
    G_k = cell(nmax,1);
    F_k = cell(nmax,1);
    D = Polyhedron('A',Gc,'b',Fc);
    D.minHRep();
    for i = 1:nmax
        G_k{i} = D.A;
        F_k{i} = D.b;
    end
end

%% Controlled Invariant Set in Two Moves
[cisLiftedA,cisLiftedb] = legacyCIS(Ac,Bc,Gc,Fc,G_k,F_k,L,method,verbose);

% Use MPT3 and let it decide.
P = Polyhedron('A',cisLiftedA,'b',cisLiftedb);
P.minHRep;
Px = P.projection(1:n,'ifourier');  
% 'ifourier' seems to be better than 'mplp' for many cases.

% If there are input constraints, project from the extended space to the
% Brunovsky space.
if (isExtended)
    Px = Px.projection(1:(n-m));
end
cisA = Px.A;
cisb = Px.b;
% Return to original coordinates:
cisA = cisA*Pmat;

disp('Projection done!')