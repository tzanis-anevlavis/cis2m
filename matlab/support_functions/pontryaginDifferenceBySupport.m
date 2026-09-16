function [G_difference, F_difference] = pontryaginDifferenceBySupport(G, F, disturbance_map, W)
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
    % Compute the Pontryagin difference between an H-polyhedron
    %       S = {x | G x <= F}
    % and a linear image of a nonempty bounded MPT Polyhedron W:
    %       S_difference = S - disturbance_map * W.
    %
    % For each row g_i of G,
    %       S_difference = {x | g_i*x <= F_i - support_W((g_i * disturbance_map)')}.
    %
    % The inequality normals remain unchanged. The formula is exact even when
    % S_difference is lower dimensional or empty, and when W is lower dimensional.

    validateattributes(G, {'numeric'}, {'2d', 'real', 'finite'}, ...
        mfilename, 'G', 1);
    validateattributes(F, {'numeric'}, {'2d', 'real', 'finite'}, ...
        mfilename, 'F', 2);
    validateattributes(disturbance_map, {'numeric'}, ...
        {'2d', 'real', 'finite'}, mfilename, 'disturbance_map', 3);

    if (size(F, 2) ~= 1 || size(F, 1) ~= size(G, 1))
        error('cis2m:pontryaginDifferenceBySupport:InvalidBounds', ...
            'F must be a column vector with one entry per row of G.');
    end
    if (size(disturbance_map, 1) ~= size(G, 2))
        error('cis2m:pontryaginDifferenceBySupport:InvalidMapRows', ...
            'disturbance_map must have one row per column of G.');
    end
    if (~isa(W, 'Polyhedron') || ~isscalar(W))
        error('cis2m:pontryaginDifferenceBySupport:InvalidDisturbance', ...
            'W must be a scalar MPT Polyhedron.');
    end
    if (W.isEmptySet() || ~W.isBounded())
        error('cis2m:pontryaginDifferenceBySupport:InvalidDisturbance', ...
            'W must be nonempty and bounded.');
    end
    if (size(disturbance_map, 2) ~= W.Dim)
        error('cis2m:pontryaginDifferenceBySupport:InvalidMapColumns', ...
            'disturbance_map must have one column per dimension of W.');
    end

    G_difference = G;
    F_difference = F;
    directions = G * disturbance_map;
    active_rows = any(directions ~= 0, 2);
    if (any(active_rows))
        support = W.support(full(directions(active_rows, :)'));
        F_difference(active_rows) = F(active_rows) - support;
    end
end
