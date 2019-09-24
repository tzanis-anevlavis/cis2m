%% 3D Example of "Computing controlled invariant sets in two moves"
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
%
%
%
%% Description:
% This script generates an example in \R^3, in which there no constraints
% on the input u. 
%
% Original system: x+ = Ao x + Bo
% State constraints: Go x <= Fo
%
% This function makes use of the Multi-Parametric Toolbox 3.0:
% M. Herceg, M. Kvasnica, C. Jones, and M. Morari, 
% ``Multi-Parametric Toolbox 3.0,'' in Proc. of the European Control 
% Conference, Zürich, Switzerland, July 17-19 2013, pp. 502-510, 
% http://control.ee.ethz.ch/ mpt.

%% Exammple Setup:
% Original system and polyhedron:
Ao = [0 1 -2; 3 -4 5; -6 7 8]; Bo = [-1; 2; 4];

Go = [-0.9497    0.3212   -0.8258;
     -0.1578   -0.2323    0.6042;
     -0.6318    0.2547    0.9783;
	  0.4516   -0.9567   -0.8661;
     -0.2593    0.8211    0.8788;
      0.6831    0.6011   -0.9636;
      0.4685    0.4917    0.3677;
      0.1421    0.6262    0.5675;
     -0.6463   -0.2334    0.0683;
      0.9148    0.2346    0.7707;
     -0.4694    0.1510    0.7980;
      0.8492    0.0601    0.2519;
     -0.5525   -0.4499   -0.7243;
     -0.2529   -0.5027   -0.5644;
     -0.8250   -0.0967   -0.6357;
      0.2802   -0.5446   -0.9164;
     -0.6388    0.6089   -0.7861;
     -0.9099    0.9722    0.2329;
      0.4463   -0.9400    0.8793;
     -0.3051    0.0713   -0.2911];

Fo = [  0.4106;
        0.9843;
        0.9456;
        0.6766;
        0.9883;
        0.7668;
        0.3367;
        0.6624;
        0.2442;
        0.2955;
        0.6802;
        0.5278;
        0.4116;
        0.6026;
        0.7505;
        0.5835;
        0.5518;
        0.5836;
        0.5118;
        0.0826];
 
D = Polyhedron('H',[Go Fo]);

A = Ao; B = Bo;
G = Go; F = Fo;

n = size(A,1);
m = size(G,1);
% Bring system in Brunovsky normal form space:
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
Ac = [zeros(n-1,1) eye(n-1); zeros(1,n)];
Bc = [zeros(n-1,1); 1];
disp('..done!')

% Compute controlled invariant set in two moves:
[mcisA,mcisb] = mcisCF(Ac,Bc,Gc,F);
[cisA,cisb,guard] = jProject(mcisA,mcisb,n);
cisMat = [cisA*Pmat cisb];
cis = Polyhedron('H', cisMat);
cis = cis.minHRep;
% Compute MCIS using MPT3:
% (note: depending on the polyhedron D it might not converge)
system = LTISystem('A',Ao,'B',Bo);
mcisEx = system.invariantSet('X',D,'maxIterations',100);
% Plotting:
figure; plot(D, 'color', 'blue', mcisEx, 'color', 'lightgray', cis, 'color', 'white')
figure; plot(mcisEx, 'color', 'lightgray', cis, 'color', 'white')
figure; plot(cis,'color','white')

D1 = D.projection(2:3);
cis1 = cis.projection(2:3);
mcisEx1 = mcisEx.projection(2:3);
figure; plot(D1, 'color', 'blue',mcisEx1, 'color', 'lightgray', cis1, 'color', 'white')
D2 = D.projection([1,3]);
cis2 = cis.projection([1,3]);
mcisEx2 = mcisEx.projection([1,3]);
figure; plot(D2, 'color', 'blue', mcisEx2, 'color', 'lightgray', cis2, 'color', 'white')
D3 = D.projection(1:2);
cis3 = cis.projection(1:2);
mcisEx3 = mcisEx.projection(1:2);
figure; plot(D3, 'color', 'blue', mcisEx3, 'color', 'lightgray', cis3, 'color', 'white')
