function [A,B,E,G,F] = inputExt(A,Pmat,Ac,Bc,Ec,Gc,F,Gu,Fu)
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
% Input constraints:
% If there are input constraints, we extend the system by one dimension to
% incorporate them.

n = size(Ac,1);	% dimension of the system
m = size(Bc,2);	% number of inputs

% Matrix A in Brunovsky form before feedback:
tmpA = (Pmat*A)/Pmat;
if (m==1)           % Single input case.
    alpha = tmpA(end,:);
    alpha_e = [-alpha 1];   % -a^T x + v = [-a^T 1] [x,v]
    % Extended system:
    Ae = [Ac Bc; zeros(1,size(Ac,2)+size(Bc,2))];
    Be = [zeros(size(Ac,1),1); 1];
    Ee = [Ec;0];
    % Extended safe set:
    Ge = [Gc zeros(size(Gc,1),size(Gu,2)); Gu.*alpha_e];
    Fe = [F; Fu];
    A = Ae; B = Be; E = Ee; G = Ge; F = Fe;
    % Size of extended system:
    n = n + 1;
    
elseif (m>1)        % Multi-input case.
    error('Coming soon');
end