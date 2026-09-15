function [G_k, F_k] = construct_Sk(Ac, Bc, Gc, Fc, Ec, Gw, Fw, q, nmax)
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
%% Description:
%
% For a safe set S and a disturbance set W, this function:
%   a)  Computes the accumulated disturbance set W_t of a linear system for 
%       monotonically increasing `t` up to and including the Minimal Robust 
%       Positively Invariant Subset (RPIS).
%   b)  Computes the Pontryagin (Minkowski) differences representing shrunk safe sets:
%                   S_t = S - W_t, t >= 1, with S_0 = S,
%       where W_t = \sum_{i=1}^t Ac^{i-1} Ec W, and the minimal RPIS
%       W_{\infty} = \sum_{i=1}^nmax Ac^{i-1} Ec W.
%   c)  If W is empty, S_t = S for all t.
%
% Inputs:   Ac, Bc, Ec : matrices that define the discrete-time linear system:
%                           x+ = Ac x + Bc u + Ec w,
%                        in the Brunovsky normal form.
%           Gc, Fc: matrices that define the safe set: 
%                           S = {x \in \R^n | Gc x <= Fc}.
%           Gw, Fw: matrices that define the disturbance set: 
%                           W = {w \in \R^k | Gw w <= Fw}.
%                   If no disturbance use: Ec = [], Gw = [], and Fw = [].
%           q:      a positive integer, q = tau + lambda, the total sequence length
%           nmax:   the largest controllability index of the system.
%
% Outputs:  G_k, F_k: (nmax + q)-by-1 cell arrays of inequality matrices
%                   and right-hand sides, respectively. 
%                   For each actual time t = 0, ..., nmax + q - 1, 
%                   cell t + 1 represents:
%                       G_k{t + 1} * [x; v] <= F_k{t + 1},
%                   i.e., S_t x R^(m*q), where n = size(Ac, 2),
%                   m = size(Bc, 2), and v contains the virtual inputs.
%                   With r_t inequalities, G_k{t + 1} is r_t-by-(n + m*q)
%                   and F_k{t + 1} is r_t-by-1. The last m*q columns of
%                   G_k{t + 1} are zero, leaving v unconstrained here.
%                   Cell 1 represents S_0 = S; cell nmax + 1 represents
%                   S_nmax = S - W_infty. All subsequent cells repeat it.
%                   If Ec is empty, every cell represents S x R^(m*q).
%
% This script makes use of the Multi-Parametric Toolbox 3.0:
% M. Herceg, M. Kvasnica, C. Jones, and M. Morari,
% ``Multi-Parametric Toolbox 3.0,'' in Proc. of the European Control
% Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510,
% http://control.ee.ethz.ch/mpt.

n = size(Ac,2); 
m = size(Bc,2);

% Basic sanity checks (fail fast if assumptions are violated)
if (~isscalar(nmax) || nmax < 1 || floor(nmax) ~= nmax)
    error('nmax must be a positive integer.');
end
if (~isscalar(q) || q < 1 || floor(q) ~= q)
    error('q must be a positive integer.');
end
if (size(Ac, 1) ~= n)
    error('Ac must be square.');
end
if (size(Bc, 1) ~= n)
    error('Rows of Bc and Ac do not match.');
end
if (size(Gc, 2) ~= n)
    error('Columns of Gc and Ac do not match.');
end
if (~isempty(Ec) && size(Ec, 1) ~= n)
    error('Rows of Ec and Ac do not match.');
end

N = nmax + q;
mq = m * q;

% Initialize collection of sets {G_k, F_k}, k = 0 .. (N - 1).
G_k = cell(N, 1);
F_k = cell(N, 1);

% Precompute base constraints. This is safe set S.
G_base = [Gc sparse(size(Gc, 1), mq)];
F_base = Fc;
% S_0 = S. Note: since MATLAB is 1-indexed, notice that 
% t = 0 corresponds to index 1. 
G_k{1} = G_base;
F_k{1} = F_base;

% Optional performance knob for complex disturbance iterations.
% Disabled by default to preserve baseline behavior.
enable_min_hrep = false;
min_hrep_every = 5;
min_hrep_min_constraints = 500;

if (~isempty(Ec))
    % In presence of disturbance, compute the Pontryagin difference:
    %       S_t = S - W_t = S - \sum_{i=1}^t Ac^{i-1} Ec W, for t >= 1
    S = Polyhedron('A', Gc, 'b', Fc);
    W = Polyhedron('A', Gw, 'b', Fw);
    A_curr = speye(n);

    % First stage: build up to S - W_{\infty} at index (nmax + 1)
    % instead of index nmax since MATLAB is 1-indexed.
    last_shrink_idx = nmax + 1;
    for t = 1:(last_shrink_idx - 1)
        % Compute S_t
        S = S - full(A_curr * Ec) * W; % MPT needs the full version.
        if (enable_min_hrep && (mod(t, min_hrep_every) == 0) && (size(S.A, 1) >= min_hrep_min_constraints))
            S.minHRep();
        end
        G_k{t + 1} = [S.A sparse(size(S.A, 1), mq)];
        F_k{t + 1} = S.b;
        A_curr = A_curr * Ac;
    end

    % Second stage: keep S - W_{\infty} for the remaining indices
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
