function [mcisA,mcisb] = mcisCF(Ac,Bc,Gc,F)

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
% This function takes as input a system in Brunovsky normal form 
% and a polyhedron in the corresponding coordinates, and returns the
% Maximal Controlled Invariant Set (MCIS) (in the higher dimensional space)
%
% System:                   x+ = Ac x + Bc u
% Polyhedral constraint:    D = {x \in \R | Gc x <= F}
% 
% Returns: matrices mcisA, mcisb such that
%                   MCIS = {x| mcisA x <= mcisb}
%
%
% This function makes use of the Multi-Parametric Toolbox 3.0:
% M. Herceg, M. Kvasnica, C. Jones, and M. Morari, 
% ?Multi-Parametric Toolbox 3.0,? in Proc. of the European Control 
% Conference, Zürich, Switzerland, July 17?19 2013, pp. 502? 510, 
% http://control.ee.ethz.ch/ mpt.

n = size(Ac,1);
m = size(Gc,1);

%% Construct MCIS in High Dimension
disp('Computing MCIS in closed-form, in high dimensional space..')
% Construct Ghat0:
Ghat0 = zeros(n*m,n+n*m);
Hhat0 = zeros(m,n+n*m);
% "Index" matrices:
P = zeros(n*m,1);
T = zeros(n*m,1);
for j = 1:m
    for i = 1:n
        % Ghat0:
        Ghat0((j-1)*n+i,i) = Gc(j,i);
        Ghat0((j-1)*n+i,j*n+i) = -1;
        Hhat0(j,j*n+i) = 1;
        % "Index" mats:
        if (Gc(j,i) > 0)
            P(n*(j-1)+i) = -1/Gc(j,i);
        elseif (Gc(j,i) < 0)
            T(n*(j-1)+i) = 1/Gc(j,i);
        end
    end
end
% Construct the combinatorial constraints matrix
Gamma = [];
Pidx = find(P~=0);
Tidx = find(T~=0);
if (n>5)    % use parallelization if we are in dimension > 5.
    parfor p = 1:size(Pidx,1)
        for t = 1:size(Tidx,1)
            temp = zeros(1,n+m*n);
            temp(n+Pidx(p)) = P(Pidx(p));
            temp(n+Tidx(t)) = T(Tidx(t));
            Gamma = [Gamma; temp];
        end
    end
else        % else it is fast enough.
    for p = 1:size(Pidx,1)
        for t = 1:size(Tidx,1)
            temp = zeros(1,n+m*n);
            temp(n+Pidx(p)) = P(Pidx(p));
            temp(n+Tidx(t)) = T(Tidx(t));
            Gamma = [Gamma; temp];
        end
    end
end
    
% Costruct the iterative constraints matrix Gtilde
Ahat = [Ac zeros(n,n*m); zeros(n*m,n) eye(n*m,n*m)];
Bhat = [Bc; zeros(n*m,1)];
Gtilde = [];
% Ridx1 = indices that involve 'u'.
% tempA1 = Ghat0 Ahat^{1} till Ghat0 Ahat^{n} at each step 'i' without
% indices in Ridx1.
% tempB1 = Ghat0 Ac^{0} Bc till Ghat0 Ac^{n-1} Bc at each step 'i' without
% indices in Ridx1.
% Gtilde = concatinate vertically the above resulting inequalities.
tempA1 = Ghat0*Ahat;
tempB1 = Ghat0*Bhat;
for i = 1:n
    Ridx1 = (tempB1~=0);
    tempA1(Ridx1,:) = [];
    Gtilde = [Gtilde; tempA1];
    tempB1 = tempA1*Bhat;
    tempA1 = tempA1*Ahat;
end
Ftilde = zeros(size(Gtilde,1),1);

% MCIS in closed-form:
mcisHmat = [Hhat0 F; Ghat0 zeros(size(Ghat0,1),1); Gtilde Ftilde; Gamma zeros(size(Gamma,1),1)];
disp('..done!')

P = Polyhedron('H',mcisHmat);
% disp(size(P.H))

% Remove redundant inequalities (there will be many depending on the dim)
disp('Obtaining minimum representation..')
P = P.minHRep;
% disp(size(P.H))
disp('..done!')

mcisHmat = P.H;
% Return matrices mcisA, mcisb: MCIS = {x| mcisA c <= mcisb}
mcisA = mcisHmat(:,1:end-1);
mcisb = mcisHmat(:,end);