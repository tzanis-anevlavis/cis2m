function tests = test_validate_brunovsky_normal_form
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    matlab_dir = fileparts(fileparts(mfilename('fullpath')));
    testCase.addTeardown(@path, path);
    addpath(fullfile(matlab_dir, 'support_functions'));
end

function testAcceptsCanonicalChainsAndNumericalResiduals(testCase)
    % Exercise unequal multi-input chains in both dense and sparse storage.
    chain_lengths = [4 2 1];
    [Ac, Bc] = brunovskyChains(chain_lengths);
    [actual_lengths, actual_nmax] = validateBrunovskyNormalForm(Ac, Bc);
    verifyEqual(testCase, actual_lengths, chain_lengths);
    verifyEqual(testCase, actual_nmax, 4);

    [actual_lengths, actual_nmax] = validateBrunovskyNormalForm( ...
        sparse(Ac), sparse(Bc));
    verifyEqual(testCase, actual_lengths, chain_lengths);
    verifyEqual(testCase, actual_nmax, 4);

    % Small factorization residuals around canonical zeros and ones are valid.
    Ac_perturbed = Ac;
    Ac_perturbed(1, 1) = 1e-11;
    Bc_perturbed = Bc;
    Bc_perturbed(4, 1) = 1 + 1e-11;
    [actual_lengths, actual_nmax] = validateBrunovskyNormalForm( ...
        Ac_perturbed, Bc_perturbed);
    verifyEqual(testCase, actual_lengths, chain_lengths);
    verifyEqual(testCase, actual_nmax, 4);

    % A scalar fully actuated system is the smallest valid canonical pair.
    [actual_lengths, actual_nmax] = validateBrunovskyNormalForm(0, 1);
    verifyEqual(testCase, actual_lengths, 1);
    verifyEqual(testCase, actual_nmax, 1);
end

function testRejectsInvalidDimensions(testCase)
    verifyError(testCase, @() validateBrunovskyNormalForm(zeros(2, 3), ...
        zeros(2, 1)), 'cis2m:validateBrunovskyNormalForm:NonSquareAc');
    verifyError(testCase, @() validateBrunovskyNormalForm(zeros(2), ...
        zeros(3, 1)), 'cis2m:validateBrunovskyNormalForm:InvalidBcRows');
    verifyError(testCase, @() validateBrunovskyNormalForm(zeros(2), ...
        zeros(2, 3)), 'cis2m:validateBrunovskyNormalForm:TooManyInputs');
end

function testRejectsMalformedInputMatrix(testCase)
    [Ac, Bc] = brunovskyChains([3 2]);
    invalid_matrices = {
        replaceEntry(Bc, 3, 2), ...       % A chain endpoint must equal one.
        'cis2m:validateBrunovskyNormalForm:InvalidInputMatrix';
        replaceEntry(Bc, 1, 1e-6), ...    % Off-endpoint entries must be zero.
        'cis2m:validateBrunovskyNormalForm:InvalidInputMatrix';
        replaceEntry(Bc, 3, -1), ...      % Endpoint signs cannot be negative.
        'cis2m:validateBrunovskyNormalForm:InvalidInputMatrix';
        zeros(size(Bc)), ...              % Every input must terminate a chain.
        'cis2m:validateBrunovskyNormalForm:InvalidInputMatrix';
        [0 0; 1 1; 0 0; 0 0; 0 0], ...  % Chain endpoints must be distinct.
        'cis2m:validateBrunovskyNormalForm:InvalidChainOrdering';
        [0 0; 0 1; 1 0; 0 0; 0 0], ...  % Endpoints must follow input order.
        'cis2m:validateBrunovskyNormalForm:InvalidChainOrdering';
        [1 0; 0 1; 0 0; 0 0; 0 0], ...  % The chains must cover every state.
        'cis2m:validateBrunovskyNormalForm:InvalidChainOrdering'
    };

    for index = 1:size(invalid_matrices, 1)
        verifyError(testCase, ...
            @() validateBrunovskyNormalForm(Ac, invalid_matrices{index, 1}), ...
            invalid_matrices{index, 2});
    end
end

function testRejectsMalformedStateMatrix(testCase)
    [Ac, Bc] = brunovskyChains([3 2]);
    invalid_matrices = {
        replaceEntry(Ac, sub2ind(size(Ac), 1, 2), 0);  % A required shift is missing.
        replaceEntry(Ac, sub2ind(size(Ac), 1, 2), 2);  % Shift gains must equal one.
        replaceEntry(Ac, sub2ind(size(Ac), 2, 2), 1);  % Diagonal dynamics are forbidden.
        replaceEntry(Ac, sub2ind(size(Ac), 3, 4), 1);  % Separate chains cannot be coupled.
        replaceEntry(Ac, sub2ind(size(Ac), 5, 1), 1);  % Terminal rows must be zero.
        replaceEntry(Ac, sub2ind(size(Ac), 1, 1), 1e-6) % Residual exceeds tolerance.
    };

    for index = 1:numel(invalid_matrices)
        verifyError(testCase, ...
            @() validateBrunovskyNormalForm(invalid_matrices{index}, Bc), ...
            'cis2m:validateBrunovskyNormalForm:InvalidStateMatrix');
    end
end

function matrix = replaceEntry(matrix, index, value)
    matrix(index) = value;
end

function [Ac, Bc] = brunovskyChains(chain_lengths)
    n = sum(chain_lengths);
    m = numel(chain_lengths);
    Ac = zeros(n);
    Bc = zeros(n, m);
    first = 1;
    for channel = 1:m
        last = first + chain_lengths(channel) - 1;
        for state = first:(last - 1)
            Ac(state, state + 1) = 1;
        end
        Bc(last, channel) = 1;
        first = last + 1;
    end
end
