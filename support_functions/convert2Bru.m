function [Ac,Bc,Ec,Gc,Fc,Pmat,Mu,nmax,alpha,Bm,guard] = convert2Bru(A,B,E,G,F,Gu,Fu)
%% Authors: Tzanis Anevlavis.
% Copyright (C) 2019, Tzanis Anevlavis.
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
% For any comments contact Tzanis Anevlavis @ janis10@ucla.edu.
%
%
%
%
%% Description:
% Convert system in Brunovsky normal form space:

n = size(A,1);  % dimension of the system
m = size(B,2);  % number of inputs

% Controllability matrix:
C = B;
for i = 1:n-1
    C = [B A*C];
end
if (rank(C)<n)
    warning('System not controllable');
    disp('System dimension:');
    disp(n);
    disp('Controllability rank:');
    disp(rank(C));
end

% Controllability Indices
Mu = zeros(m,1);
count = 0;
for i = 1:n
    for j = 1:m
        count = count+1;
        if (count==1)
            S = C(:,1); 
            Mu(j)= Mu(j)+1;
        else
            if (rank([S C(:,count)]) > rank(S))
                S = [S C(:,count)]; Mu(j)= Mu(j)+1;
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
Cbar = [];
for v = 1:m
    for i = 1:Mu(v)
        Cbar = [Cbar A^(i-1)*B(:,v)];
    end
end
% Check condition number prior to using Cbar^{-1}:
if (cond(Cbar)>1e14)
    warning('Condition number > 1e14')
end

% Similarity transformation matrix:
CbarInv = eye(n)/Cbar;  % more stable than inv(Cbar).
Pmat = [];
q = [];
for v = 1:m
    q = [q;CbarInv(sigma(v),:)];
    for i = 1:Mu(v)
        Pmat = [Pmat; q(v,:)*A^(i-1)];
    end
end

% Bm
Bm = [];
for i = 1:m
    Bm = [Bm; q(i,:)*A^(Mu(i)-1)];
end
Bm = Bm*B;
% alpha
alpha = [];
tmpA = (Pmat * A) / Pmat;
for i = 1:m
    alpha = [alpha; tmpA(sigma(i),:)];
end
% u = -inv(Bm)*ALPHA x + inv(Bm) v

% System in Brunovsky Normal Form after feedback:
Ac = (Pmat * A) / Pmat - ((Pmat*B)/Bm)*alpha;
Bc = (Pmat*B)/Bm;
Ec = Pmat*E;
nmax = max(Mu);

% Domain in Brunovsky coordinates:
Gc = G/Pmat;
Fc = F;

%% Input constraints:
% If there are input constraints, we extend the system by one dimension to
% incorporate them.

guard = false;
if (~isempty(Gu))
    guard = true;
    % u = -inv(Bm)alpha z + inv(Bm)v = inv(Bm)[-alpha I] [z,v]
    alpha_e = [-(eye(m)/Bm)*alpha eye(m)/Bm];
    % Extended system:
    Ae = [Ac Bc; zeros(size(Bc,2),size(Ac,2)+size(Bc,2))];
    Be = [zeros(size(Ac,1),size(Bc,2)); eye(size(Bc,2))];
    Ee = [Ec; zeros(size(Bc,2),size(Ec,2))];
    % Extended safe set:
    Ge = [Gc zeros(size(Gc,1),size(Gu,2)); Gu*alpha_e];
    Fe = [Fc; Fu];
    
    Ac = Ae; Bc = Be; Ec = Ee; Gc = Ge; Fc = Fe;
    
    nmax = nmax+1;
end
