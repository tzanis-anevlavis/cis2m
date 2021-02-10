function [cisA, cisb] = computeCIS(A,B,Gx,Fx,T,L,Gu,Fu,E,Gw,Fw,method,verbose)
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
%           method: if not specified the default (best for most scenarios)
%           algorithm will be selected -- legacy options for reference to
%           our papers:
%                       CDC19, HSCC20, ACC21a, ACC21b.
%
%              L : denotes the level of the hierarchy (L-invariance).
%              T : transient before L (TL-invariance).
%
%           verbose = 0 - no messages; 1 - displays messages.
%
% Output:   [cisA, cib]: matrix defining the computed Controlled Invariant Set
%                   s.t. CIS = {x \in \R^n | cisA x <= cisb}

%% Add support folder to path
% addpath('./support_functions/');

%% Input arguments check
inputArgCheck

%% Convert system in Brunovsky normal form space, and extend state space:
[Ac,Bc,Ec,Gc,Fc,Pmat,~,nmax,~,~,isExtended] = convert2Bru(A,B,E,Gx,Fx,Gu,Fu);
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
if (strcmp(method,'default'))
    [cisLiftedA,cisLiftedb] = closedformCIS(Ac,Bc,Gc,Fc,G_k,F_k,T,L,nmax,verbose);
else
    [cisLiftedA,cisLiftedb] = legacyCIS(Ac,Bc,Gc,Fc,G_k,F_k,L,method,verbose);
end

% Use MPT3 and let it decide.
P = Polyhedron('A',cisLiftedA,'b',cisLiftedb);
P.minHRep;
Px = P.projection(1:n,'ifourier');  
% 'ifourier' seems to be better than 'mplp' for many cases.

%% Comment out for debugging:
% if (isExtended)
%     % Check if result is empty:
%     if (Px.isEmptySet)
%         warning('result is empty');
%     end
%     % Check if result is contained by the safe set:
%     Sc = Polyhedron('H',[Gc Fc]);
%     if (~(Px <= Sc))
%         warning('Out of safe set.');
%     end
%     % Check if result is numerically invariant in Brunovsky extended space:
%     if (~isInvariant(Px,[],[],Ac,Bc))
%         warning('Result not numerically invariant.');
%     else
%         disp('Invariant in Brunovsky (extended) space!')
%     end
% end

%%
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

%% Comment out for debugging:
% % Check if result is empty:
% if (Px.isEmptySet)
%     warning('result is empty');
% end
% % Check if result is contained by the safe set:
% Gbru = Gx/Pmat;
% Fbru = Fx;
% Dbru = Polyhedron('H',[Gbru Fbru]);
% if (~(Px <= Dbru))
%     warning('Out of safe set');
% end
% % Check if result is numerically invariant:
% [~,~,~,~,~,~,~,~,alpha,Bm,~] = convert2Bru(A,B,E,Gx,Fx,Gu,Fu);
% Abru = (Pmat * A) / Pmat - ((Pmat*B)/Bm)*alpha;
% Bbru = (Pmat*B)/Bm;
% XU = Polyhedron('H', [Gc(size(Gbru,1)+1:end,:) Fc(size(Gbru,1)+1:end)]);
% if (~isInvariant(Px,[],XU,Abru,Bbru))
%     warning('Result not numerically invariant');
% else
%     disp('Invariant in Brunovsky space!')
% end