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
% Input arguments check.

% If no verbose, choose silent mode.
if (~exist('verbose','var'))
    verbose = 0;
end

% Check state constraints and system:
if (~exist('Gx','var') || ~exist('Fx','var'))
    error('State constraints not specified.');
end
if (size(Gx,1)~=size(Fx,1))
    error('Rows of Gx and Fx do not match.')
elseif (size(Gx,2)~=size(A,2))
    error('Columns of A (number of states) and Gx do not match.')
elseif (size(A,1)~=size(B,1))
    error('Rows of A and B do not match.')
end

% Check input constraints:
if (~exist('Gu','var') && ~exist('Fu','var'))
    Gu = [];    Fu = [];
elseif (exist('Gu','var') && ~exist('Fu','var'))
    error('Matrix Gu given, but not vector Fu.')
elseif (~exist('Gu','var') && exist('Fu','var'))
    error('Vector Fu given, but not matrix Gu.')
elseif (~isempty(Gu) && ~isempty(Fu))
    if (size(Gu,1)~=size(Fu,1))
        error('Rows of Gu and Fu do not match.')
    elseif (size(B,2)~=size(Gu,2))
        error('Columns of B (number of inputs) and Gu do not match.')
    end
end

% Disturbance matrix:
if (~exist('E','var') || isempty(E))
    disturbance = false;
    E = zeros(size(A,2),1); Gw = []; Fw = [];
else
    disturbance = true;
    if (~exist('Gw','var') && ~exist('Fw','var'))
        error('Disturbance matrix is specified, but not disturbance set.')
    elseif (exist('Gw','var') && ~exist('Fw','var'))
        error('Matrix Gw given, but not vector Fw.')
    elseif (~exist('Gw','var') && exist('Fu','var'))
        error('Vector Fw given, but not matrix Gw.')
    elseif (~isempty(Gw) && ~isempty(Fw))
        if (size(Gw,1)~=size(Fw,1))
            error('Rows of Gw and Fw do not match.')
        elseif (size(E,2)~=size(Gw,2))
            error('Columns of E and Gw do not match.')
        else
            W = Polyhedron('H',[Gw Fw]);
            if (~W.isBounded)
                error('Disturbance is unbounded.')
            end
        end
    end
end

% If no method selected, use the default:
if (~exist('method','var'))
    method = 'default'; 
end

% Method compatibility with disturbance case:
if(disturbance)
    if ( (strcmp(method,'CDC19')) || (strcmp(method,'HSCC20')) || (strcmp(method,'ACC21a')) )
        error('The selected method has not been implemented for the case of disturbances. Please select ACC21b.')
    end
end

% Loop selection:
if (~(exist('L','var') && (L>0) && (floor(L)==L) ) )
    disp('Length of the loop (L) must be a positive integer.')
    prompt = 'Choose a value for length of loop (L>0):';
    L = input(prompt);
    while ((L<0) || (L==0) || (floor(L)~=L))
        prompt = 'Input L must be a positive integer:';
        L = input(prompt);
    end
end

if (strcmp(method,'CDC19'))
    L = 1;
end