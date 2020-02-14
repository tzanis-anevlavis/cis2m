function [cisMat] = CIS2M_hierarchy(A,B,G,F,Gu,Fu,L,verbose)

%% Authors: Tzanis Anevlavis, Paulo Tabuada
% Copyright (C) 2019, Tzanis Anevlavis, Paulo Tabuada
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
% This code is part of the implementation of the algorithm proposed in:
% Tzanis Anevlavis and Paulo Tabuada, "A simple hierarchy for computing
% controlled invariant sets", Submitted to 23rd ACM International 
% Conference on Hybrid Systems: Computation and Control (HSCC'20),
% and is publicly available at: https://github.com/janis10/cis2m-hierarchy
%
% For any comments contact Tzanis Anevlavis @ janis10@ucla.edu.
%
%
%
%
%% Description:
% This function takes as input a discete-time linear system and a
% polyhedral safe set and computes a Controlled Invariant Subset of the
% safe set.
%
% Inputs:   A, B : matrices defining the discrete-time LTI: x+ = Ax + Bu
%           G, F : matrices defining the safe set: {x \in \R^n | G x <= F}
%           Gu,Fu: matrices defining the input constr: {u \in \R | Gu x <= Fu}
%                  if there are no input constraints use Gu = [], Fu = [].
%           L    : positive integer defining the length of the FTS loop
%
%           verbose = 0 - no messages; 1 - displays messages.
%
% Output:   cisMat: matrix defining the computed Controlled Invariant Set
%                   cisMat = [cisA cib] s.t. CIS = {x \in \R^n | cisA x <= cisb}

%% Add support folder to path
addpath('./support_functions/');

%% Input arguments check
if (~exist('verbose','var'))
    verbose = 0;
end
if (~exist('Gu','var') && ~exist('Fu','var'))
    Gu = [];    Fu = [];
elseif (exist('Gu','var') && ~exist('Fu','var'))
    error('Vector Fu given, but not matrix Gu.')
elseif (~exist('Gu','var') && exist('Fu','var'))
    error('Matrix Gu given, but not vector Fu.')
end
if (nargin<4)
    error('Not enough input arguments.')
end
if (~exist('L','var'))
    warning('Variable L not specified, assuming L = 1.')
    L = 1;
end


n = size(A,1);	% dimension of the original system
m = size(B,2);	% dimension of the input
k = size(G,1);  % number of constraints

%% Convert system in Brunovsky normal form space:
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
Gc(abs(Gc)<max(abs(Gc))/1e12) = 0;
% Make all elements rational numbers with 7-decimal-digit precision.
Gc = Gc.*1e7;
Gc = round(Gc);
Gc = Gc./1e7;
% System in Brunovsky Normal Form after feedback:
if (m==1)   % Single-input case.
    Ac = [zeros(n-1,1) eye(n-1); zeros(1,n)];
    Bc = [zeros(n-1,1); 1];
else        % Multi-input case.
    error('Coming soon');
end

%% Input constraints:
% If there are input constraints, we extend the system by one dimension to
% incorporate them.
guard = false;
if (~isempty(Gu))
    guard = true;
    
    % Matrix A in Brunovsky form before feedback:
    tmpA = (Pmat*A)/Pmat;
    if (m==1)           % Single input case.
        alpha = tmpA(end,:);
        alpha_e = [-alpha 1];   % -a^T x + v = [-a^T 1] [x,v]
        % Extended system:
        Ae = [Ac Bc; zeros(1,size(Ac,2)+size(Bc,2))];
        Be = [zeros(size(Ac,1),1); 1];
        % Extended safe set:
        Ge = [Gc zeros(size(Gc,1),size(Gu,2)); Gu.*alpha_e];
        Fe = [F; Fu];
        
        A = Ae; B = Be; G = Ge; F = Fe;
        % Size of extended system:
        n = n + 1;
        
    elseif (m>1)        % Multi-input case.
        error('Coming soon');
    end
else
    A = Ac; B = Bc; G = Gc;
end

%% Controlled Invariant Set in Two Moves - Hierarchy
if (L==1)
    % If L=1 - use the optimized code from CDC2019 algorithm.
    [cisLiftedA,cisLiftedb] = mcisCF(A,B,G,F,verbose);
else
    % else use our extended hierarchy code:
	[cisLiftedA,cisLiftedb] = mcisFTS(A,B,G,F,L,verbose);
end
[cisA,cisb] = jProject(cisLiftedA,cisLiftedb,n,verbose);
cisMat = [cisA cisb];

% If we had input constraints, we project from the extended space to the
% original space, i.e., eliminate input u.
if (guard)
    cisMat = fourier(cisMat,1:n-1);
end

% Return to original coordinates:
cisMat = [cisMat(:,1:end-1)*Pmat cisMat(:,end)];

%% Remove support folder from path
rmpath('./support_functions/');