function [Gstate,Ginput,Gvirtual,flifted] = computeImplicitCIS(A,B,Gx,Fx,L,Gu,Fu,E,Gw,Fw,method,verbose)
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
% This is the wrapper of the CIS2M repository.
%
% Takes as input a discete-time linear system and a polyhedral safe set.
% Computes a  (robust) Controlled Invariant Subset of the safe set.
%
% Inputs:   A, B : matrices defining the discrete-time LTI: x+ = Ax + Bu
%           Gx,Fx: matrices defining the safe set: {x \in \R^n | G x <= F}
%           Gu,Fu: matrices defining the input constr: {u \in \R | Gu u <= Fu}
%                  if there are no input constraints use Gu = [], Fu = [].
%              E : disturbance matrix if the system is: x+ = Ax + Bu + Ew
%           Gw,Fw: matrices defining the disturb. set: {w \in \R^k | Gw w <= Fw}
%                  if there is no disturbance use E = [], Gw = [], Fw = [].
%
%              L : lentgh of the loop for HSCC20, ACC21a, ACC21b (L>0 int.)
%
%           verbose = 0 - no messages; 1 - displays messages.
%
% Output:   [Gstate,Ginput,Gvirtual,flifted] : representation of the
%                   implicit controlled invariant set in the form:
%                   Gstate x + Ginput u + Gvirtual v <= flifted. 
%                   with x \in \R^n, u \in \R^m, and v \in \R^{L x m}. 

%% Add support folder to path
% addpath('./support_functions/');

%% Input arguments check
inputArgCheck

%% Convert system in Brunovsky normal form space, and extend state space:
[Ac,Bc,Ec,Gc,Fc,Pmat,~,nmax,alpha,Bm,isExtended] = convert2Bru(A,B,E,Gx,Fx,Gu,Fu);
A = Ac; B = Bc; E = Ec; Gx = Gc; Fx = Fc;
n = size(A,2);
m = size(B,2);

%% Construct D_k
% This is for ACC21b only.
if (disturbance)
    % Minkowski difference for the closed-form expression.
    G_k = cell(nmax,1);
    F_k = cell(nmax,1);
    W = Polyhedron('A',Gw,'b',Fw);
    D = Polyhedron('A',Gx,'b',Fx);
    for i = 1:nmax
        D = D - A^(i-1)*E*W;
        D.minHRep();
        G_k{i} = D.A;
        F_k{i} = D.b;
    end
else
    G_k = cell(nmax,1);
    F_k = cell(nmax,1);
    D = Polyhedron('A',Gx,'b',Fx);
    D.minHRep();
    for i = 1:nmax
        G_k{i} = D.A;
        F_k{i} = D.b;
    end
end

%% Controlled Invariant Set in One Move
if (strcmp(method,'default'))
    [cisLiftedA,cisLiftedb] = closedformCIS(Ac,Bc,Gc,Fc,G_k,F_k,L,nmax,verbose);
else
    [cisLiftedA,cisLiftedb] = legacyCIS(Ac,Bc,Gc,Fc,G_k,F_k,L,method,verbose);
end

if (isExtended)
    % Extract state input matrices in [z,u,v], z\in\R^n, u\in\R^m,
    % v\in\R^mL
    n = n-m;

    Gz = cisLiftedA(:,1:n);
    Gv = cisLiftedA(:,n+1:n+m);
    Gvirtual = cisLiftedA(:,n+m+1:end);
    % up to this point it checks out as invariant wrt:
    %       mcis = Polyhedron('H',[cisLiftedA cisLiftedb];
    %       Gmcis = mcis.A;
    %       fmcis = mcis.b;
    %       Pre = Polyhedron('H', [Gmcis*A_hd fmcis]);
    %       mcis <= Pre
    
    % Map back from Brunovsky to original space
    Gstate = Gz*Pmat+Gv*alpha*Pmat;
    Ginput = Gv*Bm;
    
    % ok and the set:
%     mcis2 = Polyhedron('H',[Gstate Ginput Gvirtual Fcl]);
%     is invariant wrt:
%     A_hd_OG = [ A                      B                    zeros(n,L*m); 
%                -inv(Bm)*alpha*Pmat*A  -inv(Bm)*alpha*Pmat*B T; 
%                            zeros(L*m,n+m)                   P]
%     
%     Gmcis = mcis2.A;
%     fmcis = mcis2.b;
%     Pre = Polyhedron('H', [Gmcis*A_hd_OG fmcis]);
%     mcis2 <= Pre
else
    % Extract state input matrices in y = (z,v)
    Gz = cisLiftedA(:,1:n);
    Gvirtual = cisLiftedA(n+1,end);
    % Map back from Brunovsky to original space
    Gstate = Gz*Pmat;
    Ginput = Gvirtual(:,1:m);
    Gvirtual = Gvirtual(:,m+1:end);
end    
flifted = cisLiftedb;     
