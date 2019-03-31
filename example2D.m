%% 2D Example of "Computing controlled invariant sets in two moves"
%
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
% This script generates an example in \R^2, in which there also constraints
% on the input u. To incorporate these constraints we extend the original 
% state space by one dimension, obtaining the state y = (x,u), and then 
% introduce a new unconstrained input v governing the evolution of the 
% state u according to u_{k+1} = ?_{k}.
%
% Original system: x+ = Ao x + Bo
% State constraints: Go x <= Fo
% Input constraints: u \in [-2,2].
%
% This function makes use of the Multi-Parametric Toolbox 3.0:
% M. Herceg, M. Kvasnica, C. Jones, and M. Morari, 
% ?Multi-Parametric Toolbox 3.0,? in Proc. of the European Control 
% Conference, Zürich, Switzerland, July 17?19 2013, pp. 502? 510, 
% http://control.ee.ethz.ch/ mpt.

%% Exammple Setup:
% Original system and polyhedron:
Ao = [1 0.1; 0 1]; Bo = [0; 0.05];
Go =[  -0.7238   -0.5117;
        0.2558   -0.8091;
       -0.0528    0.7725;
       -0.1060   -0.7190;
       -0.1252    0.1868;
        0.7232   -0.9371;
        0.4235    0.6708];
Fo = [  0.2990;
        0.0983;
        0.0276;
        0.1202;
        0.0348;
        0.0921;
        0.0240];
Dsys = Polyhedron('H',[Go Fo]);

% Extended system with input:
A = [1 0.1 0; 0 1 0.05; 0 0 0]; B = [0; 0; 1];

% Extended polyhedron:
ext =[ -0.7238   -0.5117    0    0.2990;
        0.2558   -0.8091    0    0.0983;
       -0.0528    0.7725    0    0.0276;
       -0.1060   -0.7190    0    0.1202;
       -0.1252    0.1868    0    0.0348;
        0.7232   -0.9371    0    0.0921;
        0.4235    0.6708    0    0.0240;
        0         0         1    2;
        0         0        -1    2];

Dext = Polyhedron('H',ext);
Dext.minHRep;
G = Dext.A;
F = Dext.b;

n = size(A,1);
m = size(G,1);

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

%% Compute Controllled Invariant Sets:
% Compute controlled invariant set in two moves:
[mcisA,mcisb] = mcisCF(Ac,Bc,Gc,F);
[cisA,cisb,guard] = jProject(mcisA,mcisb,n);
cisMat = [cisA*Pmat cisb];
% This is in the extended space with input u.
cisExt = Polyhedron('H', cisMat);
cisExt = cisExt.minHRep;
% Eliminate input u to obtain CIS in the original space:
cis = cisExt.projection((1:n-1),'ifourier');

% Compute MCIS using MPT3: 
% (note: depending on the polyhedron Dsys it might not converge)
system = LTISystem('A',Ao,'B',Bo);
% Polyhedron of input constraints:
U = Polyhedron('H',[1 2; -1 2]);
mcisEx = system.invariantSet('X',Dsys,'U',U);

%% Plotting:
figure; plot(Dsys, 'color', 'red', mcisEx, 'color', 'lightblue', cis, 'color', 'lightgreen')
figure; plot(cis,'color','lightgreen')