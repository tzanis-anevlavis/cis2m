function [RCIS, A_hd] = computeRCIS(A,B,Gx,Fx,Gu,Fu,E,Gw,Fw,implicit,L,T)
%% Authors: Tzanis Anevlavis
% Copyright (C) 2021, Tzanis Anevlavis
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
% For any comments contact Tzanis Anevlavis @ t.anevlavis@ucla.edu.
%
%
%
%
%% Description:
% This is the wrapper of the CIS2M repository.
%
% Takes as input a discete-time linear system and a polyhedral safe set.
% Computes a  (Robust) Controlled Invariant Subset of the safe set.
%
% Inputs:   A,B,E : matrices defining the discrete-time linear system:
%                           x+ = Ax + Bu + Ew.
%           Gx,Fx: define the safe set: 
%                           {x \in \R^n | Gx x <= Fx}.
%           Gu,Fu: define the input constraints: 
%                           {u \in \R | Gu u <= Fu}.
%                  If no input constraints, use: Gu = [] and Fu = [].
%           Gw,Fw: define the disturbance set: 
%                           {w \in \R^k | Gw w <= Fw}.
%                  If no disturbance use: E = [], Gw = [], and Fw = [].
%
%           implicit:   1 - computes high-dimensional implicit (R)CIS.
%                       0 - computes explicit (R)CIS via projection.
%
%           If only L is specified, then L : L-th level of hierarchy;
%           If both L and T are specified:
%           L (lambda): loop of eventually periodic input sequence.
%           T (tau):    transient of eventually periodic input sequence.
%
% Output:   Polyhedron object RCIS.
%           If implicit = 0, then:
%                   explicit RCIS = {x \in \R^n | rcisA x <= rcisb}.
%           If implicit = 1, then:
%                   implicit RCIS = 
%           {(x,u,v) \in \R^n x \R^m x \R^{m(T+L)} | rcisA(x,u,v) <= rcisb}.

%% Add support folder to path.
% addpath('./support_functions/');

%% Input arguments check.
inputArgCheck
% Use sparse matrices for faster computations.
A = sparse(A); B = sparse(B); Gx = sparse(Gx); Fx = sparse(Fx);
Gu = sparse(Gu); Fu = sparse(Fu);
E = sparse(E); Gw = sparse(Gw); Fw = sparse(Fw);

if (exist('T','var')) % T(tau) is specified, compute (R)CIS_(Tau,Lambda). 
    Level = T+L;
else
    Level = L;        % Compute the hierarchy level (R)CIS_(L). 
end

%% Convert system in Brunovsky normal form space and extend state space.
% TODO: extended space not needed with the latest formulation. 
% Should optimize the code. 
[Ac,Bc,Ec,Gc,Fc,Pmat,nmax,Am,Bm,isExtended] = convert2Bru(A,B,E,Gx,Fx,Gu,Fu);

%% Construct S_k sets.
[G_k, F_k] = construct_Sk(Ac,Bc,Gc,Fc,Ec,Gw,Fw,Level,nmax);

%% Implicit Controlled Invariant Set.
if (exist('T','var')) % T(tau) is specified, compute (R)CIS_(Tau,Lambda). 
    [rcisLiftedA,rcisLiftedb, A_hd, K, P] = implicitclosedformRCIS(Ac,Bc,G_k,F_k,L,T,nmax);
    implicitRCIS = Polyhedron('H',[rcisLiftedA,rcisLiftedb]);
else                  % Compute the hierarchy level (R)CIS_(L). 
    for Lambda = 1:Level
        Tau = Level-Lambda;
        [rcisLiftedA,rcisLiftedb, A_hd, K, P] = implicitclosedformRCIS(Ac,Bc,G_k,F_k,Lambda,Tau,nmax);
        implicitRCIS(Lambda) = Polyhedron('H',[rcisLiftedA,rcisLiftedb]);
    end
end

%% Output:
for i = 1:length(implicitRCIS)
    rcisLiftedA = implicitRCIS(i).A;
    rcisLiftedb = implicitRCIS(i).b;
    
    if (implicit)       % Return implicit (R)CIS.
        if (isExtended)
            % Extract state input matrices in [z,u,v], z\in\R^n, u\in\R^m,
            % v\in\R^mL.
            n = size(Ac,2)-size(Bc,2); % Dimension of original space. 
            m = size(Bc,2);
            Gz = rcisLiftedA(:,1:n);
            Gv = rcisLiftedA(:,n+1:n+m);
            Gvirtual = rcisLiftedA(:,n+m+1:end);

            % Map back from Brunovsky to original space
            Gstate = Gz * Pmat + Gv * Am * Pmat;
            Ginput = Gv * Bm;
        else
            n = size(Ac,2); m = size(Bc,2);

            % Extract state input matrices in y = (z,v)
            Gz = rcisLiftedA(:,1:n);
            Gvirtual = rcisLiftedA(:,n+1:end);
            % Map back from Brunovsky to original space
            Gstate = Gz * Pmat;
            Ginput = Gvirtual(:,1:m);
            Gvirtual = Gvirtual(:,m+1:end);
        end
        rcisA = [Gstate Ginput Gvirtual];
        rcisb = rcisLiftedb;

        % Transform A_hd from Brunovsky space to original space.
        % x+ = Ax + Bu
        %    u = -inv(Bm) Am T x + inv(Bm) v 
        % => u+ = -inv(Bm) Am T A x -inv(Bm) Am T B u + inv(Bm) K {virtual}
        % and {virtual}+ = P virtual. 
        A_hd = [A B sparse(size(A,1), size(Gvirtual,2));
                -Bm\Am*Pmat*A -Bm\Am*Pmat*B Bm\K;
                sparse(size(P,1),size(A,2)+size(B,2)) P];

    else                % Return explicit (R)CIS. 
        % Use MPT3 to project back to the original space.
        rcisLifted = Polyhedron('A',rcisLiftedA,'b',rcisLiftedb);
        rcis = rcisLifted.projection(1:size(A,2),'ifourier');  % 'ifourier' seems to be better than 'mplp' for many cases.
        rcisA = rcis.A; rcisb = rcis.b;
        % Return to original coordinates:
        rcisA = rcisA * Pmat;

        disp('Projection done!')
    end
    
    RCIS(i) = Polyhedron('H',[rcisA,rcisb]);
end
