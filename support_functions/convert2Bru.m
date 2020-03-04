function [Ac,Bc,Ec,Gc,Pmat] = convert2Bru(A,B,E,G)
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

n = size(A,1);	% dimension of the system
m = size(B,2);	% number of inputs

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

% Compute the similarity transformation matrix:
if (m == 1)         % Single-input case.
    Cbar = C;
    initPos = 1;
else                % Multi-input case.
    initPos = (1:m);
    Cbar = B;
    for i = 1:(size(C,2)-m)
        j = mod(i,m);
        tmpC = Cbar;
        tmpC = [tmpC(:,1:(initPos(j+1)-1)) C(:,m+i) tmpC(:,initPos(j+1):end)];
        if (rank(tmpC)>rank(Cbar))
            Cbar = tmpC;
            initPos(j+1:end) = initPos(j+1:end)+1;
            if (rank(Cbar)==n)
                break;
            end
        end
    end
end
% Check condition number prior to using Cbar^{-1}:
if (cond(Cbar)>1e14)
    warning('Condition number > 1e14')
end

% Compute controllability indices and \sigmas:
contrInd = zeros(m,1);
sigma = zeros(m,1);
for v = 1:m
    if (v == m)
        contrInd(v) = n+1 - initPos(v);
    else
        contrInd(v) = initPos(v+1)-initPos(v);
    end
    
    if (v==1)
        sigma(v) = contrInd(v);
    else
        sigma(v) = sigma(v-1) + contrInd(v);
    end
end

% Similarity transformation matrix:
CbarInv = inv(Cbar);
Pmat = [];
for v = 1:m
    q = CbarInv(sigma(v),:);
    for i = 1:contrInd(v)
        Pmat = [Pmat; q*A^(i-1)];
    end
end

% Domain in Brunovsky coordinates:
Gc = G/Pmat;

% System in Brunovsky Normal Form after feedback:
if (m==1)   % Single-input case.
    Ac = [zeros(n-1,1) eye(n-1); zeros(1,n)];
    Bc = [zeros(n-1,1); 1];
    Ec = Pmat*E;
else        % Multi-input case.
    error('Coming soon');
end