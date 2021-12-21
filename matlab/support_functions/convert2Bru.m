function [Ac,Bc,Ec,Gc,Fc,Pmat,nmax,Am,Bm,isExtended] = convert2Bru(A,B,E,G,F,Gu,Fu)
%% Authors: Tzanis Anevlavis.
% Copyright (C) 2021, Tzanis Anevlavis.
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
% Convert system in Brunovsky normal form space.

n = size(A,2);  % number of states
m = size(B,2);  % number of inputs

sparse(A); sparse(B);

% Controllability matrix:
Co = ctrb(A,B);
rankCo = rank(Co);
if (rankCo<n)
    warning('System not controllable.');
    disp(['System dimension:',num2str(n),'.']);
    disp(['Controllability rank:',num2str(rankCo),'.']);
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

% Similarity transformation matrix Pmat:
CbarInv = speye(n)/Cbar;  % more stable than inv(Cbar).
Pmat = zeros(n,n);
q = zeros(m,n);
idx = 1;
for v = 1:m
    q(v,:) = CbarInv(sigma(v),:);
    A_curr = speye(n);
    for i = 1:Mu(v)
        Pmat(idx,:) = q(v,:) * A_curr;
        A_curr = A_curr * A;
        idx = idx + 1;
    end
end
Pmat = sparse(Pmat);

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
tmpA = (Pmat * A) / Pmat;
for i = 1:m
    Am(i,:) = tmpA(sigma(i),:);
end
Am = sparse(Am);
% u = -inv(Bm) * Am * x + inv(Bm) * v

% System in Brunovsky Normal Form after feedback:
Ac = (Pmat * A) / Pmat - ((Pmat * B) / Bm) * Am;
Bc = (Pmat * B) / Bm;
if (~isempty(E))
    Ec = Pmat * E;
else
    Ec = [];
end
nmax = max(Mu);

% Domain in Brunovsky coordinates:
Gc = G/Pmat;
Fc = F;

%% Input constraints:
% If there are input constraints, we extend the system by one dimension to
% incorporate them.
isExtended = false;
if (~isempty(Gu))
    isExtended = true;
    % u = -inv(Bm)Am z + inv(Bm)v = inv(Bm)[-Am I] [z,v]
    alpha_e = [-(speye(m)/Bm)*Am speye(m)/Bm];
    % Extended system:
    Ae = [Ac Bc; sparse(size(Bc,2),size(Ac,2)+size(Bc,2))];
    Be = [sparse(size(Ac,1),size(Bc,2)); speye(size(Bc,2))];
    if (~isempty(Ec))
        Ee = [Ec; sparse(size(Bc,2),size(Ec,2))];
    else
        Ee = [];
    end
    % Extended safe set:
    Ge = [Gc sparse(size(Gc,1),size(Gu,2)); Gu * alpha_e];
    Fe = [Fc; Fu];
    
    Ac = Ae; Bc = Be; Ec = Ee; Gc = Ge; Fc = Fe;
    
    nmax = nmax+1;
end
