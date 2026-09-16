function validateComputeRCISInputs(A, B, E, Gxu, Fxu, Gw, Fw)
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
% Validate the system, constraint, and disturbance dimensions accepted by
% computeRCIS. Gxu constrains [x; u]. Empty E, Gw, and Fw indicate no
% disturbance. Zero input columns in Gxu indicate state-only constraints.

validateattributes(A, {'numeric'}, {'2d', 'real', 'finite', 'nonempty'}, ...
    mfilename, 'A');
validateattributes(B, {'numeric'}, {'2d', 'real', 'finite', 'nonempty'}, ...
    mfilename, 'B');
validateattributes(Gxu, {'numeric'}, {'2d', 'real', 'finite'}, ...
    mfilename, 'Gxu');
validateattributes(Fxu, {'numeric'}, {'2d', 'real', 'finite'}, ...
    mfilename, 'Fxu');

% Check the system and joint state-input constraints.
n = size(A, 1);
if (size(A, 2) ~= n)
    error('A must be square.')
elseif (size(B, 1) ~= n)
    error('Rows of A and B do not match.')
elseif (size(Gxu, 2) ~= n + size(B, 2))
    error('cis2m:validateComputeRCISInputs:InvalidJointColumns', ...
        'Gxu must have n + m columns, ordered as [x; u].');
elseif (size(Gxu, 1) ~= size(Fxu, 1) || size(Fxu, 2) ~= 1)
    error('cis2m:validateComputeRCISInputs:InvalidJointBounds', ...
        'Fxu must be a column vector with one entry per row of Gxu.');
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
    validateattributes(E, {'numeric'}, {'2d', 'real', 'finite'}, mfilename, 'E');
    validateattributes(Gw, {'numeric'}, {'2d', 'real', 'finite'}, mfilename, 'Gw');
    validateattributes(Fw, {'numeric'}, {'2d', 'real', 'finite'}, mfilename, 'Fw');
    W = Polyhedron('A', Gw, 'b', Fw);
    if (W.isEmptySet())
        error('cis2m:validateComputeRCISInputs:EmptyDisturbance', ...
            'Disturbance set must be nonempty. Use E = [], Gw = [], Fw = [] for no disturbance.');
    elseif (~W.isBounded)
        error('Disturbance set is unbounded.')
    end
end

end
