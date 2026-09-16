function [G_k, F_k] = constructShrunkSafeSetConstraints(Ac, Bc, Gc, Fc, Ec, Gw, Fw, q, nmax)
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
% This code is part of the implementation of the algorithm proposed in:
%
% For any comments contact Tzanis Anevlavis @ t.anevlavis@ucla.edu.
%
%
%% Description:
%
% For a safe set S \in \R^n x \R^m, i.e., joint in the space of states and inputs,
% and a disturbance set W, this function:
%   a)  Computes the accumulated disturbance set W_t of a linear system for
%       monotonically increasing `t` up to and including the Minimal Robust
%       Positively Invariant Subset (RPIS).
%
%   b)  Computes the Pontryagin (Minkowski) differences representing shrunk safe sets:
%                   S_t = S - (W_t x {0}), t >= 1, with S_0 = S,
%       where W_t = \sum_{i=1}^t Ac^{i-1} Ec W, and the minimal RPIS
%       W_{\infty} = \sum_{i=1}^nmax Ac^{i-1} Ec W by nilpotency of Ac.
%       Notice that only the state component is disturbed.
%
%   c)  If W is empty, S_t = S for all t.
%
% Inputs:   Ac, Bc, Ec : matrices that define the discrete-time linear system:
%                           z+ = Ac z + Bc r + Ec w,
%                        in the Brunovsky normal form.
%           Gc, Fc: matrices that define the safe set:
%                           S = {(z, r) | Gc * [z; r] <= Fc}.
%                   Gc has n + m columns, where n = size(Ac, 2)
%                   and m = size(Bc, 2).
%           Gw, Fw: matrices that define the disturbance set:
%                           W = {w \in \R^k | Gw w <= Fw}.
%                   If no disturbance use: Ec = [], Gw = [], and Fw = [].
%           q:      a positive integer, q = tau + lambda, the total sequence length
%           nmax:   the largest controllability index of the system.
%                   It must match the chain structure of (Ac, Bc).
%
% Outputs:  G_k, F_k: (nmax + q)-by-1 cell arrays of inequality matrices
%                   and right-hand sides, respectively.
%                   For each actual time t = 0, ..., nmax + q - 1,
%                   cell t + 1 represents:
%                       G_k{t + 1} * [z; r] <= F_k{t + 1}.
%                   All G_k entries equal Gc. Only the right-hand sides
%                   change: each row loses the support of W_t along its
%                   state normal. Pure r constraints are unchanged.
%                   Cell 1 represents S_0 = S; cell nmax + 1 represents
%                   S_nmax = S - (W_infty x {0}). Later cells repeat it.
%                   These are joint constraints; r = H*v is substituted by
%                   computeImplicitClosedFormRCIS for each lasso component.
%
% This script makes use of the Multi-Parametric Toolbox 3.0:
% M. Herceg, M. Kvasnica, C. Jones, and M. Morari,
% ``Multi-Parametric Toolbox 3.0,'' in Proc. of the European Control
% Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510,
% http://control.ee.ethz.ch/mpt.

%% Validate inputs.
n = size(Ac, 2);
m = size(Bc, 2);

validateattributes(nmax, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'}, mfilename, 'nmax');
validateattributes(q, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'}, mfilename, 'q');
validateattributes(Gc, {'numeric'}, {'2d', 'real', 'finite'}, mfilename, 'Gc');
validateattributes(Fc, {'numeric'}, {'2d', 'real', 'finite'}, mfilename, 'Fc');
nmax = double(nmax);
q = double(q);
if (size(Ac, 1) ~= n)
    error('Ac must be square.');
end
if (size(Bc, 1) ~= n)
    error('Rows of Bc and Ac do not match.');
end
[~, brunovsky_nmax] = validateBrunovskyNormalForm(Ac, Bc);
if (nmax ~= brunovsky_nmax)
    error('cis2m:constructShrunkSafeSetConstraints:InvalidNilpotencyIndex', ...
        ['nmax must equal the largest Brunovsky controllability index. ' ...
        'Expected %d, received %d.'], brunovsky_nmax, nmax);
end
if (size(Gc, 2) ~= n + m)
    error('cis2m:constructShrunkSafeSetConstraints:InvalidJointColumns', ...
        'Gc must have n + m columns, ordered as [z; r].');
end
if (size(Fc, 2) ~= 1 || size(Fc, 1) ~= size(Gc, 1))
    error('cis2m:constructShrunkSafeSetConstraints:InvalidJointBounds', ...
        'Fc must be a column vector with one entry per row of Gc.');
end
if (~isempty(Ec) && size(Ec, 1) ~= n)
    error('Rows of Ec and Ac do not match.');
end
if (~isempty(Ec))
    if (isempty(Gw) || isempty(Fw))
        error('With disturbance, Ec, Gw, and Fw must all be nonempty.');
    end
    W = Polyhedron('A', Gw, 'b', Fw);
    if (W.isEmptySet() || ~W.isBounded())
        error('cis2m:constructShrunkSafeSetConstraints:InvalidDisturbance', ...
            'Disturbance set must be nonempty and bounded.');
    end
else
    if (~isempty(Gw) || ~isempty(Fw))
        error('Without disturbance, Ec, Gw, and Fw must all be empty.');
    end
end

%% Construct the shrunk safe sets iteratively.
N = nmax + q;

% Initialize collection of sets {G_k, F_k}, k = 0 .. (N - 1).
G_k = cell(N, 1);
F_k = cell(N, 1);

% Precompute base constraints. This is safe set S.
G_base = Gc;
F_base = Fc;
% S_0 = S. Note: since MATLAB is 1-indexed, notice that
% t = 0 corresponds to index 1.
G_k{1} = G_base;
F_k{1} = F_base;

if (~isempty(Ec))
    % Each disturbance contribution acts only on z. Embed its map in joint
    % [z; r] space with zero rows for the Brunovsky input coordinates.
    A_curr = speye(n);

    % First stage: construct S_1, ..., S_nmax. At the start of iteration
    % t, A_curr = Ac^(t-1). Since cell 1 stores S_0, S_nmax is stored in
    % cell nmax + 1.
    last_shrink_idx = nmax + 1;
    for t = 1:(last_shrink_idx - 1)
        disturbance_map = [A_curr*Ec; sparse(m, size(Ec, 2))];
        [G_k{t + 1}, F_k{t + 1}] = pontryaginDifferenceBySupport(G_k{t}, F_k{t}, disturbance_map, W);
        A_curr = A_curr * Ac;
    end

    % Second stage: Ac^nmax = 0, so no later disturbance contribution can
    % enlarge W_t. Repeat S - (W_{\infty} x {0}).
    for t = (nmax + 1):(N - 1)
        G_k{t + 1} = G_k{last_shrink_idx};
        F_k{t + 1} = F_k{last_shrink_idx};
    end
else
    % If no disturbance, they are all the same.
    for t = 1:(N - 1)
        G_k{t + 1} = G_base;
        F_k{t + 1} = F_base;
    end
end
