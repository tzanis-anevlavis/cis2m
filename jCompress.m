function [A,b] = jCompress(A,b)

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
% This function checks if an inequality is always satisfied given the
% remaining inequalities. If it does, it is reduntant and we remove it.

idx = zeros(size(A,1),1);   % Indices to be removed
for i = 1:size(A,1)
    % Add as an extra constraint an upper bounf of b(i)+1 to avoid
    % returning INF
    tmpA = [A; A(i,:)];
    tmpb = [b; b(i)+1];
    % Find the max value of the LHS of the inequality 'i'
    tmpA(i,:) = [];
    tmpb(i) = [];
    [~,fval] = linprog(-A(i,:),tmpA,tmpb);
    val = -fval;
    
    if (val <= b(i))
        A(i,:) = 0;     % Use zeros so that the indexing in the loop is
        b(i) = 0;       % not affected and remove after the loop is done.
        idx(i) = i;
    end
end

idx = (idx~=0);
A(idx,:) = [];
b(idx,:) = [];