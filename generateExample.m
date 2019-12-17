%% Authors: Tzanis Anevlavis, Paulo Tabuada
% Copyright (C) 2019, Tzanis Anevlavis, Paulo Tabuada
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful;
% See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
% This code is part of the implementation of the algorithm proposed in:
% Tzanis Anevlavis and Paulo Tabuada, "Computing controlled invariant sets
% in two moves", in 2019 IEEE Conference on Decision and Control,
% and is publicly available at: https://github.com/janis10/cis2m
%
% For any comments contact Tzanis Anevlavis @ janis10@ucla.edu.
%
%
%
%
%% Generate case study:
% This script generates a random, compact polyhedron in \R^n, with k
% facets, and a linear control system in Brunovsky normal form.
%
% Three required variables in workspace (A, B, m):
% Matrices A,B : LTI system in terms of matrices (A,B).
% m = number of facets for the polyhedron to be generated.
%
% This script makes use of the Multi-Parametric Toolbox 3.0:
% M. Herceg, M. Kvasnica, C. Jones, and M. Morari, 
% ``Multi-Parametric Toolbox 3.0,'' in Proc. of the European Control 
% Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510, 
% http://control.ee.ethz.ch/ mpt.

%% Generate random bounded Polyhedron - with m constraints:

% Needs dimension size 'n', and number of constraints 'k' in workspace.

disp('Generating polyhedron..')
D = Polyhedron('H',[0 1 1]);
while (~D.isBounded && ~D.isEmptySet)
    G = -1 + (1+1)*rand(k,n);
    F = rand(k,1);
    D = Polyhedron('H',[G F]);
end
% Make all elements rational numbers with 10-decimal-digit precision.
G = G.*1e10;
G = round(G);
G = G./1e10;
F = F.*1e10;
F = round(F);
F = F./1e10;
D = Polyhedron('H',[G F]);
disp('..done!')

%% and a system in Brunovsky normal form:

disp('System in Brunovsky normal form..')
C = B;
for i = 1:n-1
    C = [B A*C];
end
Cinv = inv(C);
q = Cinv(end,:);
Pmat = q;
for i = 1:n-1
    Pmat = [q; Pmat*A];
end
Pmatinv = inv(Pmat);
% Domain in Brunovsky coordinates:
Gc = G*Pmatinv;
% System in Brunovsky Normal Form:
% Ac = Pmat*A*Pmatinv;
% Bc = Pmat*B;
% % State-feedback
% Ac(end,:) = zeros(1,size(Ac,2));
Ac = [zeros(n-1,1) eye(n-1); zeros(1,n)];
Bc = [zeros(n-1,1); 1];

disp('..done!')