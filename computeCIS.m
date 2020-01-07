function [cisMat] = computeCIS(A,B,G,F,Gu,Fu,verbose)

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
% This function takes as input a discete-time linear system and a
% polyhedral safe set and computes a Controlled Invariant Subset of the
% safe set.
%
% Inputs:   A, B : matrices defining the discrete-time LTI: x+ = Ax + Bu
%           G, F : matrices defining the safe set: {x \in \R^n | G x <= F}
%           Gu,Fu: matrices defining the input constr: {u \in \R | Gu x <= Fu}
%                  if there are no input constraints use Gu = [], Fu = [].
%
%           verbose = 0 - no messages; 1 - displays messages.
%
% Output:   cisMat: matrix defining the computed Controlled Invariant Set
%                   cisMat = [cisA cib] s.t. CIS = {x \in \R^n | cisA x <= cisb}

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
elseif (size(G,1)~=size(F,1))
    error('Dimensions (number of rows) of G and F do not match.')
elseif (size(G,2)~=size(A,2))
    error('Dimensions (number of columns) of A and G do not match.')
elseif (size(A,1)~=size(B,1))
    error('Dimensions of A and B do not match.')
elseif (~isempty(Gu) && ~isempty(Fu))
    if (size(Gu,1)~=size(Fu,1))
        error('Dimensions (number of rows) of Gu and Fu do not match.')
    elseif (size(B,2)~=size(Gu,2))
        error('Dimensions (number of columns) of B and Gu do not match.')
    end
end
if (nargin<4)
    error('Not enough input arguments.')
end

n = size(A,1);	% dimension of the system

%% Convert system in Brunovsky normal form space:
C = B;
for i = 1:n-1
    C = [B A*C];
end
% Check condition number prior to using C^{-1}:
if (cond(C)>1e14)
    warning('Condition number > 1e14')
end
if (rank(C)<n)
    warning('System not controllable');
    disp('System dimension:');
    disp(n);
    disp('Controllability rank:');
    disp(rank(C));
end
Cinv = inv(C);
q = Cinv(end,:);
Pmat = q;
for i = 1:n-1
    Pmat = [q; Pmat*A];
end
% Domain in Brunovsky coordinates:
Gc = G/Pmat;
Gc(abs(Gc)<max(abs(Gc))/1e12) = 0;
% Make all elements rational numbers with 7-decimal-digit precision.
Gc = Gc.*1e7;
Gc = round(Gc);
Gc = Gc./1e7;
% System in Brunovsky Normal Form:
Ac = [zeros(n-1,1) eye(n-1); zeros(1,n)];
Bc = [zeros(n-1,1); 1];

%% Input constraints:
% If there are input constraints, we extend the system by one dimension to
% incorporate them.
guard = false;
if (~isempty(Gu))
    guard = true;
    
    % Matrix A in Brunovsky form before feedback:
    tmpA = (Pmat*A)/Pmat;
    alpha = tmpA(end,:);
    alpha_e = [-alpha 1];   % -a^T x + v = [-a^T 1] [x,v]
    
    % Extended system:
    Ae = [Ac Bc; zeros(1,size(Ac,2)+size(Bc,2))];
    Be = [zeros(size(Ac,1),1); 1];
    
    % Extended safe set:
%     Ge = [Gc zeros(size(Gc,1),size(Gu,2)); zeros(size(Gu,1),size(Gc,2)) [1;-1]];
    Ge = [Gc zeros(size(Gc,1),size(Gu,2)); Gu.*alpha_e];
    Fe = [F; Fu];
    
    A = Ae; B = Be; G = Ge; F = Fe;
    
    % Size of extended system:
    n = n + 1;
end

%% Controlled Invariant Set in Two Moves
[cisLiftedA,cisLiftedb] = mcisCF(A,B,G,F,verbose);
[cisA,cisb] = jProject(cisLiftedA,cisLiftedb,n,verbose);
cisMat = [cisA cisb];

% If we had input constraints, we project from the extended space to the
% original space, i.e., eliminate input u.
if (guard)
    cisMat = fourier(cisMat,1:n-1);
end

% Return to original coordinates:
cisMat = [cisMat(:,1:end-1)*Pmat cisMat(:,end)];