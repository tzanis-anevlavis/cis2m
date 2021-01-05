function [isIT] = isInvariant(X,U,XU,A,B)
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
% Check if a polyhedron X is invariant with respect to a linear system:
% x^+ = A x + B u, with u \in U, and (x,u) \in XU.

n = size(A,2);
m = size(B,2);
% Pure state constraints:
Gx = X.A;
Fx = X.b;
% Pure input constraints:
if (~isempty(U))
    Gu = [zeros(size(U.A,1),n) U.A];
    Fu = U.b;
else
    Gu = [];
    Fu = [];
end
% Mixed state-input constraints:    (check this one again..)
if (~isempty(XU))
    Gxu = XU.A;
    Fxu = XU.b;
else
    Gxu = [];
    Fxu = [];
end
% Concatenate constraints:
matW =  [Gx*A   Gx*B    Fx; 
         Gu             Fu;        
         Gxu            Fxu];
% Compute Pre:
W = Polyhedron('H', matW);
Pre = W.projection(1:n,'ifourier');

% % Attempt using backwards reachable set from MPT:
% system = LTISystem('A',A,'B',B);
% if (isempty(U))
%     Pre = system.reachableSet('X',X,'N',1,'direction','backward');
% else
%     Pre = system.reachableSet('X',X,'U',U,'N',1,'direction','backward');
% end

isIT = (X <= Pre);