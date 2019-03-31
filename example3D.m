%% 3D Example of "Computing controlled invariant sets in two moves"
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
% This script generates an example in \R^3, in which there no constraints
% on the input u. 
%
% Original system: x+ = Ao x + Bo
% State constraints: Go x <= Fo
%
% This function makes use of the Multi-Parametric Toolbox 3.0:
% M. Herceg, M. Kvasnica, C. Jones, and M. Morari, 
% ?Multi-Parametric Toolbox 3.0,? in Proc. of the European Control 
% Conference, Zürich, Switzerland, July 17?19 2013, pp. 502? 510, 
% http://control.ee.ethz.ch/ mpt.

%% Exammple Setup:
% Original system and polyhedron:
Ao = [0 1 -2; 3 -4 5; -6 7 8]; Bo = [-1; 2; 4];

Go =[0.8304    0.9608   -0.6607;
   -0.5279    0.7682   -0.2108;
   -0.9722   -0.5355   -0.7135;
    0.1539   -0.7305    0.8130;
    0.0238   -0.4138   -0.4418;
    0.2655    0.8882   -0.8859;
   -0.6859    0.4937   -0.0190;
   -0.1077   -0.8368    0.8729;
   -0.8887    0.0932    0.3885;
   -0.4192   -0.9132    0.3742;
   -0.1649   -0.9065    0.1010;
    0.2987    0.9829    0.9178;
    0.6361   -0.6986   -0.7053];

Fo = [  0.0366;
        0.1687;
        0.8988;
        0.0233;
        0.3298;
        0.1828;
        0.2025;
        0.1244;
        0.2920;
        0.3577;
        0.3586;
        0.3896;
        0.5930];
 
Dsys = Polyhedron('H',[Go Fo]);

A = Ao; B = Bo;
G = Go; F = Fo;

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
cis = Polyhedron('H', cisMat);
cis = cis.minHRep;

% Compute MCIS using MPT3:
% (note: depending on the polyhedron Dsys it might not converge)
system = LTISystem('A',Ao,'B',Bo);
mcisEx = system.invariantSet('X',Dsys);

%% Plotting:
figure; plot(Dsys, 'color', 'red', mcisEx, 'color', 'lightblue', cis, 'color', 'lightgreen')
figure; plot(cis,'color','lightgreen')

cis1 = cis.slice(1);
mcisEx1 = mcisEx.slice(1);
figure; plot(mcisEx1, 'color', 'lightblue', cis1, 'color', 'lightgreen')
cis2 = cis.slice(2);
mcisEx2 = mcisEx.slice(2);
figure; plot(mcisEx2, 'color', 'lightblue', cis2, 'color', 'lightgreen')
cis3 = cis.slice(3);
mcisEx3 = mcisEx.slice(3);
figure; plot(mcisEx3, 'color', 'lightblue', cis3, 'color', 'lightgreen')