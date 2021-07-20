function [G_k, F_k] = construct_Sk(Ac,Bc,Gc,Fc,Ec,Gw,Fw,L,nmax)
%% Authors: Tzanis Anevlavis
% Copyright (C) 2021, Tzanis Anevlavis
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
% For any comments contact Tzanis Anevlavis @ t.anevlavis@ucla.edu.
%
%
%
%
%% Description:
% Construct S_k sets:
% In case we have disturbances we compute the Minkowski difference:
% S_k = SafeSet - A^(k-1) * E * DisturbanceSet, k = 1,..,nmax.
% If there are no disturbances:
% S_k = SafeSet, k = 1,..,nmax.
%
% This script makes use of the Multi-Parametric Toolbox 3.0:
% M. Herceg, M. Kvasnica, C. Jones, and M. Morari,
% ``Multi-Parametric Toolbox 3.0,'' in Proc. of the European Control
% Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510,
% http://control.ee.ethz.ch/mpt.

n = size(Ac,2); m = size(Bc,2);

G_k = cell(nmax+L,1);
F_k = cell(nmax+L,1);
if (~isempty(Ec))   % In presence of disturbance:
    % Minkowski difference S - A^(i-1)*E*W at each step.
    S = Polyhedron('A',Gc,'b',Fc);
    W = Polyhedron('A',Gw,'b',Fw);
    A_curr = eye(n);    % SVD does not support sparse matrices. 
    G_k{1} = [Gc zeros(size(Gc,1), m*L)];
    F_k{1} = Fc;
    for i = 2:nmax+L
        if (i<=nmax)
            S = S - A_curr * Ec * W;
%             S.minHRep(); % check if this makes things better or worse.
            G_k{i} = [S.A sparse(size(S.A,1), m*L)];
            F_k{i} = S.b;
            A_curr = A_curr * Ac; 
        else
            G_k{i} = G_k{nmax};
            F_k{i} = F_k{nmax};
        end
    end
else                % If no disturbance, they are all the same.
    for i = 1:nmax+L
        G_k{i} = [Gc sparse(size(Gc,1), m*L)];
        F_k{i} = Fc;
    end
end