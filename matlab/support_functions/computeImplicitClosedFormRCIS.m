function [rcisLiftedA, rcisLiftedb, A_lifted, H, P] = computeImplicitClosedFormRCIS(Ac, Bc, G_k, F_k, lambda, tau, nmax)
%% Authors: T.Anevlavis, Z.Liu, N.Ozay, and P.Tabuada
% Copyright (C) 2021, T.Anevlavis, Z.Liu, N.Ozay, and P.Tabuada
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
% Computes the closed-form implicit RCIS associated with a
% (tau, lambda)-lasso input sequence for the Brunovsky-form system:
%                       z+ = Ac z + Bc r.
%
% Inputs:   Ac, Bc: matrices defining the system in Brunovsky normal form.
%           G_k, F_k: cell arrays containing the inequality matrices and
%                     right-hand sides of the shrunk safe sets. For
%                     t = 0, ..., nmax + q - 1, where q = tau + lambda,
%                     cell t + 1 represents:
%                         G_k{t + 1} * [z; r] <= F_k{t + 1}.
%                     Each matrix has n + m columns. Disturbance shrinking
%                     is already included in these joint constraints.
%           lambda: positive integer period of the lasso sequence.
%           tau: nonnegative integer transient of the lasso sequence.
%           nmax: nilpotency index of Ac. It must match the largest
%                 controllability index encoded by (Ac, Bc).
%
% Outputs:  rcisLiftedA, rcisLiftedb: inequality representation of the
%                     implicit RCIS in the lifted coordinates [z; v].
%           A_lifted: state-transition matrix of the lifted system.
%           H: projection from v to the Brunovsky-form input r.
%           P: state-transition matrix generating the lasso sequence.
%           v has m*q entries grouped by input channel, each with q samples.
%           The current input is r = H*v; the generator evolves as v+ = P*v.

%% Validate inputs
validateattributes(Ac, {'numeric'}, {'2d', 'real', 'finite'}, ...
    mfilename, 'Ac', 1);
validateattributes(Bc, {'numeric'}, {'2d', 'real', 'finite'}, ...
    mfilename, 'Bc', 2);
validateattributes(lambda, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'}, ...
    mfilename, 'lambda', 5);
validateattributes(tau, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'nonnegative'}, ...
    mfilename, 'tau', 6);
validateattributes(nmax, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'}, ...
    mfilename, 'nmax', 7);

n = size(Ac, 2);
m = size(Bc, 2);
lambda = double(lambda);
tau = double(tau);
nmax = double(nmax);
q = tau + lambda;
N = nmax + q;

if (size(Ac, 1) ~= n)
    error('cis2m:computeImplicitClosedFormRCIS:NonSquareAc', ...
        'Ac must be square.');
end
if (size(Bc, 1) ~= n)
    error('cis2m:computeImplicitClosedFormRCIS:InvalidBcRows', ...
        'Ac and Bc must have the same number of rows.');
end
[~, brunovsky_nmax] = validateBrunovskyNormalForm(Ac, Bc);
if (nmax ~= brunovsky_nmax)
    error('cis2m:computeImplicitClosedFormRCIS:InvalidNilpotencyIndex', ...
        ['nmax must equal the largest Brunovsky controllability index. ' ...
        'Expected %d, received %d.'], brunovsky_nmax, nmax);
end
if (~iscell(G_k) || ~iscell(F_k))
    error('cis2m:computeImplicitClosedFormRCIS:InvalidConstraintCollections', ...
        'G_k and F_k must be cell arrays.');
end
if (numel(G_k) ~= N || numel(F_k) ~= N)
    error('cis2m:computeImplicitClosedFormRCIS:InvalidConstraintCount', ...
        'G_k and F_k must each contain nmax + tau + lambda cells.');
end

original_space_input_dimension = n + m;
for t = 1:N
    if (~isnumeric(G_k{t}) || ~ismatrix(G_k{t}) || ...
            ~isreal(G_k{t}) || any(~isfinite(nonzeros(G_k{t}))) || ...
            size(G_k{t}, 2) ~= original_space_input_dimension)
        error('cis2m:computeImplicitClosedFormRCIS:InvalidConstraintMatrix', ...
            ['G_k{%d} must be a finite real numeric matrix with ' ...
             'n + m columns, ordered as [z; r].'], t);
    end
    if (~isnumeric(F_k{t}) || ~ismatrix(F_k{t}) || ...
            ~isreal(F_k{t}) || any(~isfinite(nonzeros(F_k{t}))) || ...
            size(F_k{t}, 2) ~= 1 || ...
            size(F_k{t}, 1) ~= size(G_k{t}, 1))
        error('cis2m:computeImplicitClosedFormRCIS:InvalidConstraintVector', ...
            ['F_k{%d} must be a finite real numeric column vector with ' ...
             'one entry per row of G_k{%d}.'], t, t);
    end
end

%% Construct the high-dimensional dynamical system:
% Recall, q = tau + lambda.
H_bar = [1 sparse(1, q - 1)];
P_bar = [sparse(q - 1, 1) speye(q - 1);
         sparse(1, tau) 1 sparse(1, lambda - 1)];
H = kron(speye(m), H_bar);
P = kron(speye(m), P_bar);
% The lifted dynamical system is in \R^n x \R^{m q}.
A_lifted = [Ac Bc*H;
            sparse(m * q, n) P];

%% Construct high-dimensional invariant set:
% At time t:
%   a) substitute r_t = H * v_t in the joint constraint set; and
%   b) propagate [z_t; v_t] = A_lifted^t * [z_0; v_0] back to the initial lifted state.

% MATLAB cell t + 1 therefore contains:
%       [Gz_t, Gr_t * H] * A_lifted^t,  t = 0, ..., N - 1.
constraint_blocks = cell(N, 1);
A_power = speye(size(A_lifted));
for t = 1:N
    G_lifted = [G_k{t}(:, 1:n), G_k{t}(:, (n + 1):end) * H];
    constraint_blocks{t} = G_lifted * A_power;
    A_power = A_power * A_lifted;
end

rcisLiftedA = cat(1, constraint_blocks{:});
rcisLiftedb = cat(1, F_k{:});

end
