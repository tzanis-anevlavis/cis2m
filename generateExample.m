% This code is part of the implementation of the algorithm proposed in:
% Tzanis Anevlavis and Paulo Tabuada, "Computing controlled invariant sets
% in two moves", Submitted to IEEE 58th International Conference on
% Decision and Control 2019 (CDC '19),
% and is publicly available at: https://github.com/janis10/cis2m
%
% For any comments contact Tzanis Anevlavis @ janis10@ucla.edu.
%
%
%
%% Generate case study:
% This script generates a random, compact polyhedron in \R^n, with m
% facets. After that it brings the given system to Brunovsky normal form.
%
% Three required variables in workspace (A, B, m):
% Matrices A,B : LTI system in terms of matrices (A,B).
% m = number of facets for the polyhedron to be generated.
%
% This script makes use of the Multi-Parametric Toolbox 3.0:
% M. Herceg, M. Kvasnica, C. Jones, and M. Morari, 
% ?Multi-Parametric Toolbox 3.0,? in Proc. of the European Control 
% Conference, Zürich, Switzerland, July 17?19 2013, pp. 502? 510, 
% http://control.ee.ethz.ch/ mpt.

%% Generate random bounded Polyhedron - with m constraints:

% Size of system:
n = size(A,1);

disp('Generating polyhedron..')
Dsys = Polyhedron('H',[0 1 1]);
while (~Dsys.isBounded && ~Dsys.isEmptySet)
    G = -1 + (1+1)*rand(m,n);
    F = rand(m,1);
    Dsys = Polyhedron('H',[G F]);
end
disp('..done!')

%% Bring system in Brunovsky normal form space:

disp('Bringing system to Brunovsky form..')
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
Ac = Pmat*A*Pmatinv;
Bc = Pmat*B;
% State-feedback
Ac(end,:) = zeros(1,size(Ac,2));
disp('..done!')