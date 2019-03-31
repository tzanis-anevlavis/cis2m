function [slct,growth] = jSelect(A,n)

% This code is part of the implementation of the algorithm proposed in:
% Tzanis Anevlavis and Paulo Tabuada, "Computing controlled invariant sets
% in two moves", Submitted to IEEE 58th International Conference on
% Decision and Control 2019 (CDC '19),
% and is publicly available at: https://github.com/janis10/cis2m
%
% For any comments contact Tzanis Anevlavis @ janis10@ucla.edu.
%
%
%
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


