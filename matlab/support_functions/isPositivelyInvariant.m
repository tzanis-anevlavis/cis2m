function [isInv] = isPositivelyInvariant(X,A)
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
% Check if a polyhedron X is positively invariant with respect to a linear system:
% x^+ = A x.

% Pull both inequality and equality constraints back through x+ = A*x.
Pre = Polyhedron( ...
    'A', X.A*A, 'b', X.b, ...
    'Ae', X.Ae*A, 'be', X.be);

isInv = (X <= Pre);
