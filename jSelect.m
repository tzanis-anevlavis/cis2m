function [slct,growth] = jSelect(A,n)

%% Authors: Tzanis Anevlavis, Paulo Tabuada
% Copyright (C) 2019, Tzanis Anevlavis, Paulo Tabuada
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful;
% See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
% This code is part of the implementation of the algorithm proposed in:
% Tzanis Anevlavis and Paulo Tabuada, "Computing controlled invariant sets
% in two moves", in 2019 IEEE Conference on Decision and Control, 
% and is publicly available at: https://github.com/janis10/cis2m
%
% For any comments contact Tzanis Anevlavis @ janis10@ucla.edu.
%
%
%
%
%% Description:
% This function finds the variable that will yield the best growth if
% eliminated by Fourier-Motzkin Elimination.
%
% n = number of variables we do not want to eliminate

N = size(A,1);
m = size(A,2);

growth = N^2;

for i = 1:(m-n)
    pos = sum(A(:,n+i)>0);
    neg = sum(A(:,n+i)<0);
    
    g = pos*neg - (pos+neg);
    
    if (g <= growth)
        growth = g;
        slct = n+i;
    end
end


