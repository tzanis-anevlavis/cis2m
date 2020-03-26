function [mcisA,mcisb] = cdc20a(Ac,Gc,F,L,verbose)
%% Authors: T.Anevlavis, Z.Liu, N.Ozay, and P.Tabuada (alphabetically)
% Copyright (C) 2020, T.Anevlavis, Z.Liu, N.Ozay, and P.Tabuada
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

if (verbose)
    disp('Lifting problem to compute controlled invariant set in closed-form . . .')
end

%% New lift -- corresponds to 0-step:
n = size(Ac,1);
k = size(Gc,1);

% l < x < u
H = [-eye(n); eye(n)];
R = [-eye(n) zeros(n); zeros(n) eye(n)];

% B \subseteq D
GP = (Gc>0);
GN = (Gc<0);
Gl = [Gc.*GN Gc.*GP];


G0 = [H -R; zeros(k,n) Gl];
G0 = [G0 zeros(2*n+k,L)];
F0 = [zeros(2*n,1); F];

%% Construct 1:(L-1)-step matrices:

tmpG = [];
F1L1 = [];

for s = 1:L-1
    sbar = min(s,n);
    
    Phi = [Ac^s zeros(n); zeros(n) Ac^s];
    Y = [zeros(n,s-n) [zeros(n-s,sbar); eye(sbar)] zeros(n,L-s)];
    
    tmpG = [tmpG; Gl*Phi Gc*Y];
    F1L1 = [F1L1; F];
end

G1L1 = [zeros(size(tmpG,1),n) tmpG];

%% Construct L-step matrices -- corresponds to returning to the initial hyper-rectangle:
s = L;
sbar = min(s,n);
Phi = [Ac^s zeros(n); zeros(n) Ac^s];
Y = [zeros(n,s-n) [zeros(n-s,sbar); eye(sbar)] zeros(n,L-s)];

HP = (H>0);
HN = (H<0);
Hl = [H.*HN H.*HP];

GL = [zeros(2*n,n) (Hl*Phi-R) H*Y];
FL = [zeros(2*n,1)];

%% Construct the final set:
finalSet = [G0 F0; G1L1 F1L1; GL FL];

mcisA = finalSet(:,1:end-1);
mcisb = finalSet(:,end);

if (verbose)
    disp('..computed!')
end
