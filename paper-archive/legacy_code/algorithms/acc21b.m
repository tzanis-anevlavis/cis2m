function [mcisA,mcisb] = acc21b(Ac,Bc,Gc,F,G_k,F_k,L,nmax,verbose)
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

n = size(Ac,2);
m = size(Bc,2);

%% Initial constraints:
Gx = Gc;
Gv = zeros(size(Gc,1),m*L);

for t = 1:nmax+L-1
    for i = 1:L
        ABmat{i} = zeros(n,m);
    end
    for j = 1:t
        idx = mod(t-j+1-1,L)+1;
        ABmat{idx} = ABmat{idx} + Ac^(j-1)*Bc;
    end
    tbar = min(t,nmax);
    Gv = [Gv; G_k{tbar}*cat(2, ABmat{:})];
    Gx = [Gx; G_k{tbar}*Ac^t];
    F = [F; F_k{tbar}];
end

%% Construct the final set:
mcisA = [Gx Gv];
mcisb = F;

if (verbose)
    disp('..computed!')
end
