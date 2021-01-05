% function [A,B,E,G,F,Gu,Fu,Gw,Fw,method,L,verbose] = inputArgCheck(A,B,E,G,F,Gu,Fu,Gw,Fw,method,loop,verbose)
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

% If no method selected, prompt user to select:
if (~exist('method','var'))
    prompt = 'Choose method: 1 for CDC19, 2 for HSCC20, 3 for CDC20a, 4 for CDC20b.';
    method = input(prompt);
    if (method==1)
        method = 'CDC19';
    elseif (method==2)
        method = 'HSCC20';
    elseif (method==3)
        method = 'CDC20a';
    elseif (method==4)
        method = 'CDC20b';
    else
        error('Invalid choice.');
    end
end

% Method compatibility with disturbance case:
if (exist('E','var'))
    if(~isempty(E))
        if ( (strcmp(method,'CDC19')) || (strcmp(method,'HSCC20')) || (strcmp(method,'CDC20a')) )
            error('The selected method has not been implemented for the case of disturbances. Please select CDC20b.')
        end
    end
end

% Input constraints:
if (~exist('Gu','var') && ~exist('Fu','var'))
    Gu = [];    Fu = [];
elseif (exist('Gu','var') && ~exist('Fu','var'))
    error('Vector Fu given, but not matrix Gu.')
elseif (~exist('Gu','var') && exist('Fu','var'))
    error('Matrix Gu given, but not vector Fu.')
elseif (~isempty(Gu) && ~isempty(Fu))
    if (size(Gu,1)~=size(Fu,1))
        error('Dimensions (number of rows) of Gu and Fu do not match.')
    elseif (size(B,2)~=size(Gu,2))
        error('Dimensions (number of columns) of B and Gu do not match.')
    end
end
% State constraints and system:
if (size(G,1)~=size(F,1))
    error('Dimensions (number of rows) of G and F do not match.')
elseif (size(G,2)~=size(A,2))
    error('Dimensions (number of columns) of A and G do not match.')
elseif (size(A,1)~=size(B,1))
    error('Dimensions of A and B do not match.')
end
% Disturbance matrix:
if (~exist('E','var'))
    disturbance = false;
    E = zeros(n,1); Gw = []; Fw = [];
else
    disturbance = true;
    if (~exist('Gw','var') && ~exist('Fw','var'))
        error('Disturbance matrix is specified, but not disturbance set.')
    elseif (exist('Gw','var') && ~exist('Fw','var'))
        error('Vector Fw given, but not matrix Gw.')
    elseif (~exist('Gw','var') && exist('Fu','var'))
        error('Matrix Gw given, but not vector Fw.')
    elseif (~isempty(Gw) && ~isempty(Fu))
        if (size(Gw,1)~=size(Fw,1))
            error('Dimensions (number of rows) of Gw and Fw do not match.')
        elseif (size(E,2)~=size(Gw,2))
            error('Dimensions (number of columns) of E and Gw do not match.')
        else
            W = Polyhedron('H',[Gw Fw]);
            if (~W.isBounded)
                error('Disturbance is unbounded.')
            end
        end
    end
end

% Loop selection:
if (~strcmp(method,'CDC19'))
    if (exist('L','var') && (L>0) && (floor(L)==L) )
        
    else
        disp('Length of the loop (L) must be a positive integer.')
        prompt = 'Choose a value for length of loop (L>0):';
        L = input(prompt);
        while ((L<0) || (L==0) || (floor(L)~=L))
            prompt = 'Input L must be a positive integer:';
            L = input(prompt);
        end
    end
else
    L = 1;
end