function [A,b] = jProject(Asp,bsp,n,verbose)
%% Authors: Tzanis Anevlavis
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
% This function takes as input a set {x| Ax <= b} in a higher dimensional
% space in terms of the matrices A, b, and the dimension n back to which we
% want to project the aforementioned set.
%
% It returns matrices Asp, bsp that constitute the set in \R^n.
%
% Inputs:   Asp, bsp such that {x| Asp x <= bsp} in \R^m, m > n.
%           n = dimension into which the above set is projected.
%           verbose = 0 - no messages; 1 - displays messages.
%
% Outputs: A, b such that {x| A x <= b} is the projected set in \R^n.
%
% This function is a inspired by:
% A. Simon and A. King, ``Exploiting sparsity in polyhedral analysis,''
% in Proceedings of the 12th International Conference on Static Analysis,
% ser. SAS'05. Berlin, Heidelberg: Springer-Verlag, 2005, pp. 336-351.
%
% for smart selection of variables to eliminate.

%% Project back to original space using iterative FME exploiting sparsity.

% First remove redundant inequalities. This is important to have faster
% projection.
A = full(Asp);
b = full(bsp);

if (verbose)
    disp('Obtaining minimum representation..')
end

[A,b] = jCompress(A,b);
if (verbose)
    disp('..done!')
end

if (verbose)
	disp('Begin projection back to original space..')
end

% Number of constraints.
N = size(A,1);
% Number of variables to be eliminated.
V = size(A,2);
m = V-n;
% Set limit for growth.
limit = size(A,1) + ceil(0.125*N);

% Find the variable to be eliminated, and place it last.
[f_elim,~] =  jSelect(A,n);
A = A(:,[1:f_elim-1 f_elim+1:end f_elim]);

cnt = m;

% Classic Fourier-Motzkin Elimination (FME) with smart selection of vars to
% eliminate, and intermediate removal of redundant inequalities. Takes 
% advantage of parallel pool for redundancy elimination.
while (cnt>0)
    % FME:
    proj = 1;
    if (proj==0)
        tmpP = Polyhedron('H',[A b]);
        tmpP = tmpP.projection(1:(V-1));
        A = tmpP.A;
        b = tmpP.b;
    else
        tmp = fourier([A b],1:V-1);
        A = tmp(:,1:end-1);
        b = tmp(:,end);
    end

    % Update cnt and size of A.
    cnt = cnt - 1;
    V = V - 1;
    N = size(A,1);
    if (cnt>0)
        % Find the variable to be eliminated, and place it last.
        [f_elim,growth] =  jSelect(A,n);
        A = A(:,[1:f_elim-1 f_elim+1:end f_elim]);
        
        if (N+growth>limit)
            % Remove vacuously satisfied inequalities.
            % by implementing a simple LP.
            [A,b] = jCompress(A,b);
            % Update N and limit.
            N = size(A,1);
            limit = size(A,1) + ceil(0.125*N);
            % Find the variable to be eliminated, and place it last.
            [f_elim,~] =  jSelect(A,n);
            A = A(:,[1:f_elim-1 f_elim+1:end f_elim]);
        end
    end
end

if (verbose)
	disp('..done!')
end