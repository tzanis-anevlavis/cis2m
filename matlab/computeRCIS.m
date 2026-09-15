function [RCIS, A_hd] = computeRCIS(A, B, E, Gx, Fx, Gu, Fu, Gw, Fw, options)
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
%           Gx, Fx: define the safe set:
%                           {x \in \R^n | Gx x <= Fx}.
%           Gu, Fu: define the input constraints:
%                           {u \in \R^m | Gu u <= Fu}.
%                   If no input constraints, use: Gu = [] and Fu = [].
%           Gw, Fw: define the disturbance set:
%                           {w \in \R^k | Gw w <= Fw}.
%                  If no disturbance use: E = [], Gw = [], and Fw = [].
%
% Name-value options (after the nine positional matrix arguments):
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
% Example:  computeRCIS(A, B, E, Gx, Fx, Gu, Fu, Gw, Fw, lambda=3, tau=2)
%           Requires R2021a for name=value syntax; R2019b and R2020 releases
%           can use 'lambda', 3, 'tau', 2 instead.
%
% Outputs:  RCIS: Polyhedron object, or a Polyhedron array in hierarchy mode.
%           If is_implicit = 0, then:
%                   explicit RCIS = {x \in \R^n | rcisA x <= rcisb}.
%           If is_implicit = 1, then:
%                   RCIS is represented in the lifted state-input space.
%                   Its dimension depends on q and on whether input constraints
%                   require state extension. Hierarchy mode returns one
%                   Polyhedron per (tau, lambda) pair, ordered by
%                   lambda = 1, ..., q.
%           A_hd: state-transition matrix associated with the returned implicit
%                 representation. In hierarchy mode, this is the matrix for the
%                 last component (tau = 0, lambda = q).

arguments
    A
    B
    E
    Gx
    Fx
    Gu
    Fu
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

validateComputeRCISInputs(A, B, E, Gx, Fx, Gu, Fu, Gw, Fw);

% Use sparse matrices for faster computations.
A = sparse(A);
B = sparse(B);
Gx = sparse(Gx);
Fx = sparse(Fx);
Gu = sparse(Gu);
Fu = sparse(Fu);
E = sparse(E);
Gw = sparse(Gw);
Fw = sparse(Fw);

%% Convert system in Brunovsky normal form space and extend state space.
% TODO: extended space not needed with the latest formulation.
% Should optimize the code.
[Ac, Bc, Ec, Gc, Fc, Pmat, nmax, Am, Bm, isExtended] = transformToBrunovskyNormalForm(A, B, E, Gx, Fx, Gu, Fu);

%% Construct shrunk safe sets.
[G_k, F_k] = constructShrunkSafeSetConstraints(Ac, Bc, Gc, Fc, Ec, Gw, Fw, q, nmax);

%% Implicit Controlled Invariant Set.
if (~should_compute_full_hierarchy)
    % Not full hierarchy, compute (R)CIS_(tau, lambda).
    [rcisLiftedA, rcisLiftedb, A_hd, K, P] = computeImplicitClosedFormRCIS(Ac, Bc, G_k, F_k, lambda, tau, nmax);
    implicitRCIS = Polyhedron('H',[rcisLiftedA, rcisLiftedb]);
else
    % Full hierarchy computation at level q.
    for lambda = 1:q
        tau = q - lambda;
        [rcisLiftedA, rcisLiftedb, A_hd, K, P] = computeImplicitClosedFormRCIS(Ac, Bc, G_k, F_k, lambda, tau, nmax);
        implicitRCIS(lambda) = Polyhedron('H',[rcisLiftedA, rcisLiftedb]);
    end
end

%% Output:
for i = 1:length(implicitRCIS)
    rcisLiftedA = implicitRCIS(i).A;
    rcisLiftedb = implicitRCIS(i).b;
    if (is_implicit)
        % Return implicit (R)CIS.
        if (isExtended)
            % Extract state input matrices in [z,u,v], z\in\R^n, u\in\R^m, v\in\R^(m*q).
            n = size(Ac, 2) - size(Bc, 2); % Dimension of original space.
            m = size(Bc, 2);
            Gz = rcisLiftedA(:, 1:n);
            Gv = rcisLiftedA(:, (n + 1):(n + m));
            Gvirtual = rcisLiftedA(:, (n + m + 1):end);
            % Map back from Brunovsky to original space
            Gstate = Gz * Pmat + Gv * Am * Pmat;
            Ginput = Gv * Bm;
        else
            n = size(Ac, 2);
            m = size(Bc, 2);
            % Extract state input matrices in y = (z,v)
            Gz = rcisLiftedA(:, 1:n);
            Gvirtual = rcisLiftedA(:, (n + 1):end);
            % Map back from Brunovsky to original space
            Gstate = Gz * Pmat;
            Ginput = Gvirtual(:, 1:m);
            Gvirtual = Gvirtual(:, (m + 1):end);
        end
        rcisA = [Gstate Ginput Gvirtual];
        rcisb = rcisLiftedb;

        % Transform A_hd from Brunovsky space to original space.
        % x+ = Ax + Bu
        %    u = -inv(Bm) Am T x + inv(Bm) v
        % => u+ = -inv(Bm) Am T A x -inv(Bm) Am T B u + inv(Bm) K {virtual}
        % and {virtual}+ = P virtual.
        A_hd = [A B sparse(size(A, 1), size(Gvirtual, 2));
                -Bm\Am*Pmat*A -Bm\Am*Pmat*B Bm\K;
                sparse(size(P, 1), size(A, 2) + size(B, 2)) P];

    else
        % Return explicit (R)CIS.
        % Use MPT3 to project back to the original space.
        rcisLifted = Polyhedron('A', rcisLiftedA, 'b', rcisLiftedb);
        rcis = rcisLifted.projection(1:size(A, 2), 'ifourier');  % 'ifourier' seems to be better than 'mplp' for many cases.
        rcisA = rcis.A;
        rcisb = rcis.b;
        % Return to original coordinates:
        rcisA = rcisA * Pmat;
        disp('cis2m:computeRCIS: Projection done!')
    end

    RCIS(i) = Polyhedron('H', [rcisA, rcisb]);
end
