function [RCIS, A_lifted] = computeRCIS(A, B, E, Gxu, Fxu, Gw, Fw, options)
%% Authors: Tzanis Anevlavis
% Copyright (C) 2026, Tzanis Anevlavis
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
% This is the wrapper of the CIS2M repository.
%
% Takes as input a discrete-time linear system and a polyhedral safe set.
% Computes a  (Robust) Controlled Invariant Subset of the safe set.
%
% Inputs:   A, B, E : matrices defining the discrete-time linear system:
%                           x+ = Ax + Bu + Ew.
%           Gxu, Fxu: define the joint state-input safe set:
%                           {(x, u) | Gxu * [x; u] <= Fxu}.
%                   Gxu has n + m columns.
%                   Separate constraints can be combined as
%                   Gxu = blkdiag(Gx, Gu), Fxu = [Fx; Fu].
%           Gw, Fw: define the disturbance set:
%                           {w \in \R^k | Gw w <= Fw}.
%                  If no disturbance use: E = [], Gw = [], and Fw = [].
%
% Name-value options (after the seven positional matrix arguments):
%           lambda:             positive integer loop length (default: []).
%           tau:                nonnegative integer transient length (default: 0).
%           hierarchy_level:    positive integer hierarchy level (default: []).
%           is_implicit:        true or 1 for implicit output (default: true),
%                               false or 0 for explicit output via projection.
%           A nonempty hierarchy_level takes precedence over lambda and tau;
%           ignored values are not validated. Otherwise, lambda is required.
%           Empty lambda and hierarchy_level mean "not specified".
%           The sequence length q is hierarchy_level in hierarchy mode,
%           or tau + lambda in single-component mode.
%
% Example:  computeRCIS(A, B, E, Gxu, Fxu, Gw, Fw, lambda=3, tau=2)
%           Requires R2021a for name=value syntax; R2019b and R2020 releases
%           can use 'lambda', 3, 'tau', 2 instead.
%
% Outputs:  RCIS: Polyhedron object, or a Polyhedron array in hierarchy mode.
%           If is_implicit = 0, then:
%                   explicit RCIS = {x \in \R^n | rcisA x <= rcisb}.
%           If is_implicit = 1, then:
%                   RCIS is in [x; v] coordinates, of dimension n + m*q.
%                   Virtual inputs v are grouped by channel, q entries each.
%                   Given z = T*x, r = H*v, the physical input is
%                   u = Bm\(H*v - Am*T*x).
%                   Hierarchy mode returns one Polyhedron per (tau, lambda)
%                   pair, ordered by lambda = 1, ..., q.
%           A_lifted: nominal lifted dynamics in [x; v] coordinates, including
%                     when is_implicit is false. Disturbances enter via [E; 0].
%                     For a single component, this is one sparse matrix. In
%                     hierarchy mode, this is a cell array where A_lifted{i}
%                     corresponds to RCIS(i), lambda = i, and tau = q - i.

arguments
    A
    B
    E
    Gxu
    Fxu
    Gw
    Fw
    options.lambda = []
    options.tau = 0
    options.hierarchy_level = []
    options.is_implicit = true
end

%% Add support folder to path.
% addpath('./support_functions/');

%% Validate inputs.
validateattributes(options.is_implicit, {'numeric', 'logical'}, ...
    {'scalar', 'real', 'binary'}, mfilename, 'is_implicit');
is_implicit = logical(options.is_implicit);

% Select the mode before validating parameters that it may override.
should_compute_full_hierarchy = ~isempty(options.hierarchy_level);
if (should_compute_full_hierarchy)
    validateattributes(options.hierarchy_level, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'integer', 'positive'}, ...
        mfilename, 'hierarchy_level');
    q = double(options.hierarchy_level);
elseif (~isempty(options.lambda))
    validateattributes(options.lambda, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'integer', 'positive'}, mfilename, 'lambda');
    validateattributes(options.tau, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'integer', 'nonnegative'}, mfilename, 'tau');
    lambda = double(options.lambda);
    tau = double(options.tau);
    q = tau + lambda;
else
    error('cis2m:computeRCIS:MissingMode', ...
        'Specify a nonempty lambda or hierarchy_level.');
end

validateComputeRCISInputs(A, B, E, Gxu, Fxu, Gw, Fw);

% Use sparse matrices for faster computations.
A = sparse(A);
B = sparse(B);
Gxu = sparse(Gxu);
Fxu = sparse(Fxu);
E = sparse(E);
Gw = sparse(Gw);
Fw = sparse(Fw);

%% Transform the system and joint constraints into Brunovsky coordinates.
[Ac, Bc, Ec, Gc, Fc, T, nmax] = transformToBrunovskyNormalForm(A, B, E, Gxu, Fxu);

%% Construct shrunk safe sets.
[G_k, F_k] = constructShrunkSafeSetConstraints(Ac, Bc, Gc, Fc, Ec, Gw, Fw, q, nmax);

%% Compute each lasso component and return it in original state coordinates.
if (should_compute_full_hierarchy)
    % This is an array of `lambda` values for the given hierarchy level.
    periods = 1:q;
    A_lifted = cell(size(periods));
else
    % This is a scalar corresponding to a specific `(tau, lambda)` pair.
    periods = lambda;
end

n = size(A, 2);
m = size(B, 2);
T_lift = blkdiag(T, speye(m * q));
for i = 1:numel(periods)
    lambda = periods(i);
    tau = q - lambda; % Recover `tau` for the specific `(tau, lambda)` pair.

    % Compute the lifted RCIS and companion dynamical system.
    [rcisLiftedA, rcisLiftedb, A_comp] = ...
        computeImplicitClosedFormRCIS(Ac, Bc, G_k, F_k, lambda, tau, nmax);

    % Express in original-state lifted coordinates `[x; v]`:
    % [z; v] = T_lift * [x; v] = [T, 0; 0, I] * [x; v] = [T x; v]
    rcisLifted = Polyhedron('A', rcisLiftedA * T_lift, 'b', rcisLiftedb);
    if (is_implicit)
        RCIS(i) = rcisLifted;
    else
        RCIS(i) = rcisLifted.projection(1:n, 'ifourier');
    end
    % Transform the companion dynamical system to the original-state lifted coordinates.
    A_comp_original = T_lift \ (A_comp * T_lift);
    if (should_compute_full_hierarchy)
        A_lifted{i} = A_comp_original;
    else
        A_lifted = A_comp_original;
    end
end

end
