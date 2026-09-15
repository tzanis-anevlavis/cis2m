function inputArgCheck(A, B, E, Gx, Fx, Gu, Fu, Gw, Fw)
%% Authors: Tzanis Anevlavis.
% Copyright (C) 2021, Tzanis Anevlavis.
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
% Validate the system, constraint, and disturbance dimensions accepted by
% computeRCIS. Empty Gu and Fu indicate no input constraints. Empty E, Gw,
% and Fw indicate no disturbance.

% Check the system and state constraints.
n = size(A, 1);
if (size(A, 2) ~= n)
    error('A must be square.')
elseif (size(B, 1) ~= n)
    error('Rows of A and B do not match.')
elseif (size(Gx, 2) ~= n)
    error('Columns of A (number of states) and Gx do not match.')
elseif (size(Gx, 1) ~= size(Fx, 1))
    error('Rows of Gx and Fx do not match.')
elseif (~isempty(Fx) && size(Fx, 2) ~= 1)
    error('Fx must be a column vector.')
end

% Check input constraints.
if (~isempty(Gu) && isempty(Fu))
    error('Input constraints incomplete: Matrix Gu given, but not vector Fu.')
elseif (isempty(Gu) && ~isempty(Fu))
    error('Input constraints incomplete: Vector Fu given, but not matrix Gu.')
elseif (~isempty(Gu))
    if (size(Gu, 1) ~= size(Fu, 1))
        error('Rows of Gu and Fu do not match.')
    elseif (size(Gu, 2) ~= size(B, 2))
        error('Columns of B (number of inputs) and Gu do not match.')
    elseif (size(Fu, 2) ~= 1)
        error('Fu must be a column vector.')
    end
end

% Check the disturbance matrix and set.
if (isempty(E))
    if (~isempty(Gw) || ~isempty(Fw))
        error('Without disturbance, E, Gw, and Fw must all be empty.')
    end
elseif (isempty(Gw) || isempty(Fw))
    error('With disturbance, E, Gw, and Fw must all be nonempty.')
elseif (size(E, 1) ~= n)
    error('Rows of A and E do not match.')
elseif (size(Gw, 2) ~= size(E, 2))
    error('Columns of E and Gw do not match.')
elseif (size(Gw, 1) ~= size(Fw, 1))
    error('Rows of Gw and Fw do not match.')
elseif (size(Fw, 2) ~= 1)
    error('Fw must be a column vector.')
else
    W = Polyhedron('A', Gw, 'b', Fw);
    if (~W.isBounded)
        error('Disturbance set is unbounded.')
    end
end

end
