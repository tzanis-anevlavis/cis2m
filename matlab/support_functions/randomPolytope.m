function P = randomPolytope(n,k)
%% Authors: Tzanis Anevlavis
% Copyright (C) 2020, Tzanis Anevlavis
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
% Computes a random polytope of dimension n with k constraints. 

P = Polyhedron('H',[1 1]); % Initialize: one dimensional halfspace.

while (~P.isBounded || P.isEmptySet) % repeat until P = non-empty and bounded.
   F = ones(k,1); G = rand(k,n)-0.5;
   P = Polyhedron('H', [G F]);
end