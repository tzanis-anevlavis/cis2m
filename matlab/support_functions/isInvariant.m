function [isInv] = isInvariant(X,U,A,B)
%% Authors: Tzanis Anevlavis.
% Copyright (C) 2026, Tzanis Anevlavis.
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
%% Description:
% Check if a polyhedron X is invariant with respect to a linear system:
% x^+ = A x + B u, with u \in U.

n = size(A,2);

% Pull the successor constraints back through x+ = A*x + B*u.
predecessor_A = [X.A*A X.A*B];
predecessor_b = X.b;
predecessor_Ae = [X.Ae*A X.Ae*B];
predecessor_be = X.be;

% Constrain u when an admissible-input polyhedron is provided. An empty U
% denotes an unconstrained input rather than an empty admissible-input set.
if (~isempty(U))
    predecessor_A = [predecessor_A;
                     zeros(size(U.A, 1), n) U.A];
    predecessor_b = [predecessor_b; U.b];
    predecessor_Ae = [predecessor_Ae;
                      zeros(size(U.Ae, 1), n) U.Ae];
    predecessor_be = [predecessor_be; U.be];
end

predecessor_xu = Polyhedron( ...
    'A', predecessor_A, 'b', predecessor_b, ...
    'Ae', predecessor_Ae, 'be', predecessor_be);
Pre = predecessor_xu.projection(1:n, 'ifourier');

isInv = (X <= Pre);
