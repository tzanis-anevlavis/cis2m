function [Ac, Bc, Ec, Gc, Fc, T, nmax, Am, Bm] = transformToBrunovskyNormalForm(A, B, E, Gxu, Fxu)
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
% Converts a linear system and a corresponding set of constraints into the
% Brunovsky normal form space.
% For a controllable linear system:
%           x+ = A x + B u,
% perform a coordinate change and linear state feedback:
%           z = T x
%           r = Am T x + Bm u
% so that the resulting system is:
%           z+ = Ac z + Bc r ,
% where Ac and Bc are in Brunovsky canonical (companion) form.
%
% The function applies the same transformation to any state / input constraints.
% Constraints Gxu * [x; u] <= Fxu become Gc * [z; r] <= Fc.

%% Compute the required matrices for the coordinate change and the state feedback.
% The goal is to compute the similarity transformation T, and the matrices Am and Bm.

n = size(A,2);  % number of states
m = size(B,2);  % number of inputs

sparse(A);
sparse(B);

% Controllability matrix:
Co = ctrb(A,B);
rankCo = rank(Co);
if (rankCo<n)
    error('cis2m:transformToBrunovskyNormalForm:UncontrollableSystem', ...
        ['System must be controllable. State dimension is %d, but the ' ...
        'controllability matrix has rank %d.'], n, rankCo);
end

% Controllability Indices
Mu = zeros(m,1);
count = 0;
S = zeros(n,n);
idx = 1;
for i = 1:n
    for j = 1:m
        count = count+1;
        if (count==1)
%             S = Co(:,1);
            S(:,idx) = Co(:,1);
            Mu(j)= Mu(j)+1;
            rank_S = rank(S);
            idx = idx + 1;
        else
            rank_Snew = rank([S Co(:,count)]);
            if (rank_Snew > rank_S)
%                 S = [S Co(:,count)];
                S(:,idx) = Co(:,count);
                Mu(j)= Mu(j)+1;
                rank_S = rank_Snew;
                idx = idx + 1;
            end
        end
    end
end
sigma = zeros(m,1);
sigma(1) = Mu(1);
for v = 2:m
    sigma(v) = sigma(v-1) + Mu(v);
end

% Compute Cbar = [b1..A^(mu_1-1)b1,...,bm..A^(mu_m-1)bm]
Cbar = zeros(n,n);
idx = 1;
for v = 1:m
    A_curr = speye(n);
    for i = 1:Mu(v)
        Cbar(:,idx) = A_curr * B(:,v);
        A_curr = A_curr * A;
        idx = idx + 1;
    end
end
% Check condition number prior to using Cbar^{-1}:
if (cond(full(Cbar))>1e14)
    warning('Condition number > 1e14.')
end

% Similarity transformation matrix T:
CbarInv = speye(n)/Cbar;  % more stable than inv(Cbar).
T = zeros(n,n);
q = zeros(m,n);
idx = 1;
for v = 1:m
    q(v,:) = CbarInv(sigma(v),:);
    A_curr = speye(n);
    for i = 1:Mu(v)
        T(idx,:) = q(v,:) * A_curr;
        A_curr = A_curr * A;
        idx = idx + 1;
    end
end
T = sparse(T);

% Bm
Bm = zeros(m,n);
for i = 1:m
    Bm(i,:) = q(i,:) * A^(Mu(i)-1);
end
Bm = Bm * B;
Bm = sparse(Bm);

% Am
% Am = [];
Am = zeros(m,n);
tmpA = (T * A) / T;
for i = 1:m
    Am(i,:) = tmpA(sigma(i),:);
end
Am = sparse(Am);

%% Transformation in Brunovsky normal form.
% Perform:
%           z = T x
%           r = Am T x + Bm u
% so that the resulting system is:
%           z+ = Ac z + Bc r .

% System in Brunovsky normal form after feedback.
Ac = (T * A) / T - ((T * B) / Bm) * Am;
Bc = (T * B) / Bm;
if (~isempty(E))
    Ec = T * E;
else
    Ec = [];
end
nmax = max(Mu);

% Substitute x = T^{-1} z and u = Bm^{-1} (r - Am T x)
% in the joint constraints.
Gr = Gxu(:, (n + 1):end) / Bm;
Gz = Gxu(:, 1:n) / T - Gr * Am;
Gc = [Gz Gr];
Fc = Fxu;

end
