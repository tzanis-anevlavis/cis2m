function [mcisA,mcisb] = hscc20(Ac,Bc,Gc,F,L,verbose)
%% Authors: T.Anevlavis, and P.Tabuada
% Copyright (C) 2020, T.Anevlavis, and P.Tabuada
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
% Computes a lifted controlled invariant set using Algorithm 1, proposed in
% Tzanis Anevlavis and Paulo Tabuada, "A simple hierarchy for computing
% controlled invariant sets", In Proceedings of the 23rd ACM International 
% Conference on Hybrid Systems: Computation and Control (HSCC'20).

if (verbose)
    disp('Lifting problem to compute controlled invariant set in closed-form . . .')
end

n = size(Ac,1);
k = size(Gc,1);
finTrace = (1:L);
%% Lift process:
% Original lift (G0,H0):
G0 = zeros(k*n,n);
H0 = zeros(k,n+k*n);
for j = 1:k
    G0((j-1)*n+1:j*n,:) = diag(Gc(j,:));
    for i = 1:n
        H0(j,j*n+i) = 1;
    end
end
% FTS-Lift:
L = length(finTrace);           % L is the number of states of the automaton,
% and also the number of sets we are
% constructing.
M = length(unique(finTrace));   % M is the different "identities" of \lambda
% and correspond to the convex sets between
% which the state jumps in the original
% space

%% Construct H-matrices for MCIS:
Hbar = zeros(k,n+M*k*n,L);      % 3D matrix such that each (:,:,l)
% corresponds to copy l
repH = repmat(H0(:,n+1:end),1,M);
for l = 1:L
    id = finTrace(l);       % each element of finTrace()
    % corresponds to the "idendity = s"
    % of the \lambda_s
    mask = zeros(size(repH));
    mask(:,(id-1)*(k*n)+1:(id)*(k*n)) = 1;
    temp = repH.*mask;
    Hbar(:,:,l) = [zeros(k,n) temp];
end
% Final constraints on \lambda to form rectangles.
Hcell = num2cell(Hbar,[1,2]);
Hhat = vertcat(Hcell{:});
Fhat = repmat(F,L,1);


%% Construct G-constraints - Pre: multiplying by the dynamics:
Ahat = [Ac zeros(n,M*n*k); zeros(M*n*k,n) eye(M*n*k,M*n*k)];
Bhat = [Bc; zeros(M*n*k,1)];

% Matrix Gbar: Gx < \lambda_s at each copy :
Gbar = zeros(k*n,n+M*k*n,L);    % 3D matrix such that each (:,:,l)
% corresponds to copy l
repEye = repmat(-eye(k*n),1,M);
for l = 1:L
    id = finTrace(l);       % each element of finTrace()
    % corresponds to the "idendity = s"
    % of the \lambda_s
    mask = zeros(size(repEye));
    mask(:,(id-1)*(k*n)+1:(id)*(k*n)) = eye(k*n);
    temp = repEye.*mask;
    Gbar(:,:,l) = [G0 temp];
end
% Construct a 4D matrix such that:
% (:,:,i+1,l) represents (GA)_{i}x < C \lambda of set S_l
FourD = [];
for l = 1:L
    FourD = cat(4,FourD,Gbar(:,:,l));
end
for l = 1:L
    tmpG = Gbar(:,:,l);
    tempA1 = tmpG*Ahat;
    tempB1 = tmpG*Bhat;
    for i = 1:n
        Ridx = (tempB1~=0);
        tempA1(Ridx,:) = 0;
        
        FourD(:,:,i+1,l) = tempA1;
        
        tempB1 = tempA1*Bhat;
        tempA1 = tempA1*Ahat;
    end
end

% All constraints from dynamics in Ghat 3D such that:
% Ghat(:,:,l)x < lambda are all constraints involving x in the S_l part of
% the MCIS.
Ghat = [];
for l = 1:L
    tmpGhat = FourD(:,:,1,l);
    for i = 1:n
        tmpGhat = [tmpGhat; FourD(:,:,1+i,mod(l+i-1,L)+1)];
    end
    Ghat = cat(3,Ghat,tmpGhat);
end
Zhat = zeros(size(Ghat(:,:,1),1),1);

%% Construction of Gamma-constraints:
% Given the finite trace of length L finTrace = y1..yL (1xL), create a
% matrix such that tracePerm(l,:), has the part of trace starting at
% position l, until l+n (repeats at the beginning if l+n>L).
traceMat = zeros(L,n);
PerMat = [zeros(L-1,1) eye(L-1); 1 zeros(1,L-1)];
repNum = fix(n/L);
if (repNum) % if n >= L, repeat finTrace repNum times
    for l = 1:L
        tmpTrace = (PerMat^(l-1)*finTrace')';
        trace_l = [repmat(tmpTrace,1,repNum) tmpTrace(1:(n-repNum*L))];
        traceMat(l,:) = trace_l;
    end
else        % else keep the part that is of interest
    for l = 1:L
        if (l+n-1<L)
            trace_l = finTrace(l:l+n-1);
        else
            trace_l = [finTrace(l:L) finTrace(1:(n-(L-l+1)))];
        end
        traceMat(l,:) = trace_l;
    end
end

% "Index" matrices:
% tmpG = flip(Gc',2);
% tmpG = flip(Gc,2)';
% tmpG = tmpG(:);
tmpG = Gc;
pIdx = (tmpG>0);
Ptot = tmpG.*pIdx;
tIdx = (tmpG<0);
Ttot = tmpG.*tIdx;

% For each of the L sub-traces, build the \Gamma constraints by taking all
% their permutations with repetition. Store those in a 3D matrix of the
% form (:,:,l).
Gamma = [];
for l = 1:L
    permCombs = permn(traceMat(l,:),2); % each element y from traceMat(l,:)
    % corresponds to a \lambda_l
    tracePos = (1:n);
    permPos = permn(tracePos,2); % this corresponds to the "position" of y in the trace
    
    for comb = 1:size(permCombs,1)
        % particular assignment order here does not matter because we take
        % all possible combinations into account
        pos1 = permPos(comb,1);
        pos2 = permPos(comb,2);
        mask1 = zeros(size(Ptot));
        mask1(:,n+1-pos1) = 1;
        P = Ptot.*mask1;
        P = P';
        P = P(:);
        mask2 = zeros(size(Ttot));
        mask2(:,n+1-pos2) = 1;
        T = Ttot.*mask2;
        T = T';
        T = T(:);
        
        % identity of \lambda_s involved
        id1 = permCombs(comb,1);
        id2 = permCombs(comb,2);
        
        % iterate only over non-zero elements
        Pidx = find(P~=0);
        Tidx = find(T~=0);
        % bring to 1/g form
        P(Pidx) = -1./P(Pidx);
        T(Tidx) = 1./T(Tidx);
        for p = 1:size(Pidx,1)
            for t = 1:size(Tidx,1)
                temp = zeros(1,M*k*n);
                temp((id1-1)*k*n+Pidx(p)) = P(Pidx(p));
                temp((id2-1)*k*n+Tidx(t)) = T(Tidx(t));
                
                Gamma = [Gamma; temp];
            end
        end
        
    end
end
Gamma = [zeros(size(Gamma,1),n) Gamma];
Dhat = zeros(size(Gamma,1),1);

%% Construct the final set:
% Projections of all the sets in this lifted set coincide so we work only
% using one for efficiency:
finalSet = [Ghat(:,:,1) Zhat; Hhat Fhat; Gamma Dhat];

mcisA = finalSet(:,1:end-1);
mcisb = finalSet(:,end);

if (verbose)
    disp('..computed!')
end