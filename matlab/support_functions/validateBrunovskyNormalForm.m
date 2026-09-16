function [chain_lengths, nmax] = validateBrunovskyNormalForm(Ac, Bc)
    %% Authors: Tzanis Anevlavis.
    % Copyright (C) 2026, Tzanis Anevlavis.
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
    % Verify that (Ac, Bc) follows the Brunovsky normal-form convention used by
    % CIS2M. States are grouped consecutively by input channel. Each Bc column
    % is one at the final state of its chain and zero elsewhere. Ac shifts each
    % chain state to its successor with a unit superdiagonal entry and is zero
    % everywhere else.
    %
    % Outputs:  chain_lengths: controllability index of each input channel.
    %           nmax: largest controllability index and nilpotency index of Ac.

    validateattributes(Ac, {'numeric'}, ...
        {'2d', 'real', 'finite', 'nonempty'}, mfilename, 'Ac', 1);
    validateattributes(Bc, {'numeric'}, ...
        {'2d', 'real', 'finite', 'nonempty'}, mfilename, 'Bc', 2);

    n = size(Ac, 2);
    m = size(Bc, 2);
    if (size(Ac, 1) ~= n)
        error('cis2m:validateBrunovskyNormalForm:NonSquareAc', ...
            'Ac must be square.');
    end
    if (size(Bc, 1) ~= n)
        error('cis2m:validateBrunovskyNormalForm:InvalidBcRows', ...
            'Ac and Bc must have the same number of rows.');
    end
    if (m > n)
        error('cis2m:validateBrunovskyNormalForm:TooManyInputs', ...
            'Brunovsky form requires no more input channels than states.');
    end

    % transformToBrunovskyNormalForm involves numerical matrix factorizations,
    % so accept small residuals around the canonical zero and one entries.
    tolerance = 1e-8;

    % Infer one candidate chain endpoint from the largest entry in each column,
    % then verify the complete matrix against the corresponding unit vectors.
    [~, chain_ends] = max(abs(Bc), [], 1);
    chain_ends = full(chain_ends);
    Bc_expected = sparse(chain_ends, 1:m, 1, n, m);
    if (maxAbsoluteEntry(Bc - Bc_expected) > tolerance)
        error('cis2m:validateBrunovskyNormalForm:InvalidInputMatrix', ...
            ['Each Bc column must be a unit vector at the final state of its ' ...
            'Brunovsky chain.']);
    end

    % Positive, ordered chain lengths partition all n states.
    if (any(diff(chain_ends) <= 0) || chain_ends(end) ~= n)
        error('cis2m:validateBrunovskyNormalForm:InvalidChainOrdering', ...
            ['Bc chain endpoints must occur in strictly increasing rows, and ' ...
            'the final endpoint must be the last state.']);
    end
    chain_lengths = diff([0 chain_ends]);

    % Every nonterminal chain state shifts to its successor; terminal rows have
    % no nominal state transition in Brunovsky normal form.
    shift_rows = setdiff(1:(n - 1), chain_ends);
    Ac_expected = sparse(shift_rows, shift_rows + 1, 1, n, n);
    if (maxAbsoluteEntry(Ac - Ac_expected) > tolerance)
        error('cis2m:validateBrunovskyNormalForm:InvalidStateMatrix', ...
            ['Ac must contain unit superdiagonal shifts within each Brunovsky ' ...
            'chain and zeros elsewhere.']);
    end

    nmax = max(chain_lengths);

end

function value = maxAbsoluteEntry(matrix)
    entries = nonzeros(matrix);
    if (isempty(entries))
        value = 0;
    else
        value = max(abs(entries));
    end
end
