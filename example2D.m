%% 2D Example of "Computing controlled invariant sets in two moves"
%
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
%% Description:
% This script generates an example in \R^2, in which there also constraints
% on the input u. To incorporate these constraints we extend the original 
% state space by one dimension, obtaining the state y = (x,u), and then 
% introduce a new unconstrained input v governing the evolution of the 
% state u according to u_{k+1} = u_{k}.
%
% Original system: x+ = Ao x + Bo
% State constraints: Go x <= Fo
% Input constraints: u \in [-1,1].
%
% This function makes use of the Multi-Parametric Toolbox 3.0:
% M. Herceg, M. Kvasnica, C. Jones, and M. Morari, 
% ``Multi-Parametric Toolbox 3.0,'' in Proc. of the European Control 
% Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510, 
% http://control.ee.ethz.ch/ mpt.

%% Exammple Setup:
% Original system:
Ao = [1.2 1; 0 1.2]; Bo = [0.5; 0.3];
% State constraints:
Go =[   0.9147   -0.5402;
        0.2005    0.6213;
       -0.8193    0.9769;
       -0.4895   -0.8200;
        0.7171   -0.3581;
        0.8221    0.0228;
        0.3993   -0.8788];
Fo = [  0.5566;
        0.8300;
        0.7890;
        0.3178;
        0.4522;
        0.7522;
        0.1099];
D = Polyhedron('H',[Go Fo]);
% Input constraints:
Gu = [1; -1]; Fu = [1; 1];
U = Polyhedron('H',[Gu Fu]);
% Extended system with input:
A = [Ao Bo; zeros(1,size(Ao,2)+size(Bo,2))]; B = [zeros(size(Ao,1),1); 1];
% Extended polyhedron:
extMat =[ Go zeros(size(Go,1),size(Gu,2)) Fo; zeros(size(Gu,1),size(Go,2)) Gu Fu];
De = Polyhedron('H',extMat);
De.minHRep;
G = De.A;
F = De.b;

n = size(A,1);
m = size(G,1);
% Bring system in Brunovsky normal form space:
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
Ac = [zeros(n-1,1) eye(n-1); zeros(1,n)];
Bc = [zeros(n-1,1); 1];

% Compute controlled invariant set in two moves:
[mcisA,mcisb] = mcisCF(Ac,Bc,Gc,F);
[cisA,cisb,~] = jProject(mcisA,mcisb,n);
cisMat = [cisA*Pmat cisb];
% This is in the extended space with input u.
cisExt = Polyhedron('H', cisMat);
cisExt = cisExt.minHRep;
% Eliminate input u to obtain CIS in the original space:
cis = cisExt.projection((1:n-1));
% Compute MCIS using MPT3: 
% (note: depending on the polyhedron Dsys it might not converge)
system = LTISystem('A',Ao,'B',Bo);
mcisEx = system.invariantSet('X',D,'U',U,'maxIterations',300);
% Plotting:
figure; plot(D, 'color', 'red', mcisEx, 'color', 'lightblue', cis, 'color', 'lightgreen')
