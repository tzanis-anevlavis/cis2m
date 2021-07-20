function [mcisA,mcisb] = legacyCIS(Ac,Bc,Gc,Fc,G_k,F_k,L,method,verbose)
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
% Computes the MCIS in the higher dimensional space for the methods
% described in the following papers:
%
% CDC19: 
% T.Anevlavis, and P.Tabuada. 2019. Computing controlled invariant sets in 
% two moves. In 2019 IEEE 58th Conference on Decision and Control (CDC).
% 2019. 6248-6254. DOI: 10.1109/CDC40024.2019.9029610. 
%
% HSCC20:
% T.Anevlavis, and P.Tabuada. 2020. A simple hierarchy for computing 
% controlled invariant sets. In Proceedings of the 23rd ACM International 
% Conference on Hybrid Systems: Computation and Control (HSCC' 20). 
% DOI: 10.1145/3365365.3382205. 
%
% ACC21:
% T.Anevlavis, Z.Liu, N.Ozay, and P.Tabuada. 2021. An enhanced hierarchy for 
% (robust) controlled invariance. 2021 American Control Conference (ACC), 
% 2021. 
%
%

% The below methods have been implemented only for single input systems.
if (size(Bc,2)>1)
    error('The legacy methods have been implemented only for the single input case.');
end

addpath('./algorithms/');

if (strcmp(method,'CDC19'))
    [mcisA,mcisb] = cdc19(Ac,Bc,Gc,Fc,verbose);
elseif (strcmp(method,'HSCC20'))
    [mcisA,mcisb] = hscc20(Ac,Bc,Gc,Fc,L,verbose);
elseif (strcmp(method,'ACC21a'))
    [mcisA,mcisb] = acc21a(Ac,Gc,Fc,L,verbose);
elseif (strcmp(method,'ACC21b'))
    [mcisA,mcisb] = acc21b(Ac,Gc,Fc,G_k,F_k,L,verbose);
else
    error('Invalid method.');
end

removepath('./algorithms/');