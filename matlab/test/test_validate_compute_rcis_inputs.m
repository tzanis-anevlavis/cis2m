function tests = test_validate_compute_rcis_inputs
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    matlab_dir = fileparts(fileparts(mfilename('fullpath')));
    testCase.addTeardown(@path, path);
    addpath(matlab_dir, fullfile(matlab_dir, 'support_functions'));
    assumeTrue(testCase, exist('Polyhedron', 'class') == 8 && ...
        exist('mpt_init', 'file') == 2, 'Requires MPT3.');
    mpt_init;
end

function testAcceptsValidInputs(testCase)
    data = validTestData();
    verifyWarningFree(testCase, @() validateComputeRCISInputs(data{:}));

    % No disturbance is represented by three empty arguments.
    no_disturbance = replaceArguments(data, [3 6 7], {[], [], []});
    verifyWarningFree(testCase, @() validateComputeRCISInputs(no_disturbance{:}));

    % A zero disturbance map is still a specified disturbance model.
    zero_map = replaceArgument(data, 3, zeros(size(data{3})));
    verifyWarningFree(testCase, @() validateComputeRCISInputs(zero_map{:}));
end

function testRejectsInvalidA(testCase)
    data = validTestData();
    A = data{1};
    invalid_values = {
        repmat('x', size(A)), 'A';                          % Nonnumeric values are not allowed.
        zeros([size(A) 2]), 'A';                            % A must be two-dimensional.
        A + 1i, 'A';                                        % Complex values are not allowed.
        replaceEntry(A, 1, NaN), 'A';                       % NaN values are not allowed.
        replaceEntry(A, 1, Inf), 'A';                       % Infinite values are not allowed.
        [], 'A';                                            % A must be nonempty.
        zeros(size(A, 1), size(A, 2) + 1), 'square'         % A must be square.
    };
    verifySingleArgumentRejections(testCase, data, 1, invalid_values);
end

function testRejectsInvalidB(testCase)
    data = validTestData();
    B = data{2};
    invalid_values = {
        repmat('x', size(B)), 'B';                          % Nonnumeric values are not allowed.
        zeros([size(B) 2]), 'B';                            % B must be two-dimensional.
        B + 1i, 'B';                                        % Complex values are not allowed.
        replaceEntry(B, 1, NaN), 'B';                       % NaN values are not allowed.
        replaceEntry(B, 1, Inf), 'B';                       % Infinite values are not allowed.
        [], 'B';                                            % B must be nonempty.
        zeros(size(B, 1) - 1, size(B, 2)), 'Rows of A and B' % A and B must have equal row counts.
    };
    verifySingleArgumentRejections(testCase, data, 2, invalid_values);
end

function testRejectsInvalidJointConstraintMatrix(testCase)
    data = validTestData();
    Gxu = data{4};
    invalid_values = {
        repmat('x', size(Gxu)), 'Gxu';                          % Nonnumeric values are not allowed.
        zeros([size(Gxu) 2]), 'Gxu';                            % Gxu must be two-dimensional.
        Gxu + 1i, 'Gxu';                                        % Complex values are not allowed.
        replaceEntry(Gxu, 1, NaN), 'Gxu';                       % NaN values are not allowed.
        replaceEntry(Gxu, 1, Inf), 'Gxu';                       % Infinite values are not allowed.
        zeros(size(Gxu, 1), size(Gxu, 2) - 1), 'n + m columns'; % Gxu must have one column per x and u entry.
        [], 'n + m columns';                                    % Plain [] does not retain n + m columns.
        zeros(0, size(Gxu, 2)), 'one entry per row'             % Gxu and Fxu must have equal row counts.
    };
    verifySingleArgumentRejections(testCase, data, 4, invalid_values);
end

function testRejectsInvalidJointConstraintBounds(testCase)
    data = validTestData();
    Fxu = data{5};
    invalid_values = {
        repmat('x', size(Fxu)), 'Fxu';            % Nonnumeric values are not allowed.
        zeros([size(Fxu) 2]), 'Fxu';              % Fxu must be two-dimensional.
        Fxu + 1i, 'Fxu';                          % Complex values are not allowed.
        replaceEntry(Fxu, 1, NaN), 'Fxu';         % NaN values are not allowed.
        replaceEntry(Fxu, 1, Inf), 'Fxu';         % Infinite values are not allowed.
        Fxu(1:(end - 1)), 'one entry per row';    % Gxu and Fxu must have equal row counts.
        Fxu', 'column vector';                    % Fxu must be a column vector.
        [], 'one entry per row'                   % Fxu must have one bound per Gxu row.
    };
    verifySingleArgumentRejections(testCase, data, 5, invalid_values);
end

function testRejectsInvalidE(testCase)
    data = validTestData();
    E = data{3};
    invalid_values = {
        repmat('x', size(E)), 'E';                           % Nonnumeric values are not allowed.
        zeros([size(E) 2]), 'E';                             % E must be two-dimensional.
        E + 1i, 'E';                                         % Complex values are not allowed.
        replaceEntry(E, 1, NaN), 'E';                        % NaN values are not allowed.
        replaceEntry(E, 1, Inf), 'E';                        % Infinite values are not allowed.
        zeros(size(E, 1) - 1, size(E, 2)), 'Rows of A and E'; % A and E must have equal row counts.
        zeros(size(E, 1), size(E, 2) - 1), 'Columns of E and Gw'; % E columns must match the disturbance dimension.
        [], 'all be empty'                                   % E cannot be empty while Gw and Fw are present.
    };
    verifySingleArgumentRejections(testCase, data, 3, invalid_values);
end

function testRejectsInvalidDisturbanceConstraintMatrix(testCase)
    data = validTestData();
    Gw = data{6};
    invalid_values = {
        repmat('x', size(Gw)), 'Gw';                             % Nonnumeric values are not allowed.
        zeros([size(Gw) 2]), 'Gw';                               % Gw must be two-dimensional.
        Gw + 1i, 'Gw';                                           % Complex values are not allowed.
        replaceEntry(Gw, 1, NaN), 'Gw';                          % NaN values are not allowed.
        replaceEntry(Gw, 1, Inf), 'Gw';                          % Infinite values are not allowed.
        zeros(size(Gw, 1), size(Gw, 2) - 1), 'Columns of E and Gw'; % Gw columns must match the disturbance dimension.
        zeros(size(Gw, 1) - 1, size(Gw, 2)), 'Rows of Gw and Fw'; % Gw and Fw must have equal row counts.
        [], 'all be nonempty'                                    % Gw cannot be empty while E and Fw are present.
    };
    verifySingleArgumentRejections(testCase, data, 6, invalid_values);
end

function testRejectsInvalidDisturbanceConstraintBounds(testCase)
    data = validTestData();
    Fw = data{7};
    invalid_values = {
        repmat('x', size(Fw)), 'Fw';          % Nonnumeric values are not allowed.
        zeros([size(Fw) 2]), 'Fw';            % Fw must be two-dimensional.
        Fw + 1i, 'Fw';                        % Complex values are not allowed.
        replaceEntry(Fw, 1, NaN), 'Fw';       % NaN values are not allowed.
        replaceEntry(Fw, 1, Inf), 'Fw';       % Infinite values are not allowed.
        Fw(1:(end - 1)), 'Rows of Gw and Fw'; % Gw and Fw must have equal row counts.
        [Fw Fw], 'column vector';             % Fw must be a column vector.
        [], 'all be nonempty'                 % Fw cannot be empty while E and Gw are present.
    };
    verifySingleArgumentRejections(testCase, data, 7, invalid_values);
end

function testDisturbancePresenceCombinations(testCase)
    data = validTestData();

    % E is absent, but at least one disturbance-set argument is present.
    invalid_combinations = {
        [3 6 7], {[], data{6}, []};
        [3 6 7], {[], [], data{7}};
        [3 6 7], {[], data{6}, data{7}};
        % E is present, but at least one disturbance-set argument is absent.
        [6 7], {[], []};
        [6 7], {[], data{7}};
        [6 7], {data{6}, []}
    };
    for i = 1:size(invalid_combinations, 1)
        candidate = replaceArguments(data, invalid_combinations{i, 1}, ...
            invalid_combinations{i, 2});
        verifyRejected(testCase, candidate, 'all be');
    end
end

function testEmptyJointConstraintCombinations(testCase)
    data = validTestData();
    joint_dimension = size(data{1}, 1) + size(data{2}, 2);

    % Empty constraints are valid when both arrays retain their required shape.
    empty_constraints = replaceArguments(data, [4 5], ...
        {zeros(0, joint_dimension), zeros(0, 1)});
    verifyWarningFree(testCase, @() validateComputeRCISInputs(empty_constraints{:}));

    % Emptying only one side leaves an inconsistent H-representation.
    empty_matrix = replaceArgument(data, 4, zeros(0, joint_dimension));
    verifyRejected(testCase, empty_matrix, 'one entry per row');
    empty_bounds = replaceArgument(data, 5, zeros(0, 1));
    verifyRejected(testCase, empty_bounds, 'one entry per row');

    % Plain [] does not preserve the required n + m columns for Gxu.
    plain_empty = replaceArguments(data, [4 5], {[], []});
    verifyRejected(testCase, plain_empty, 'n + m columns');
end

function testDisturbanceSetGeometry(testCase)
    data = validTestData();

    % A lower-dimensional bounded disturbance set is valid. Here W is the
    % singleton {(1,-2)}.
    point = [1; -2];
    singleton = replaceArgument(data, 7, [point; -point]);
    verifyWarningFree(testCase, @() validateComputeRCISInputs(singleton{:}));

    empty_set = replaceArgument(data, 7, -ones(size(data{7})));
    verifyRejected(testCase, empty_set, 'must be nonempty');

    % Repeated bounds on w1 leave w2 unconstrained.
    unbounded_Gw = [1 0; -1 0; 2 0; -2 0];
    unbounded_set = replaceArgument(data, 6, unbounded_Gw);
    verifyRejected(testCase, unbounded_set, 'unbounded');
end

function data = validTestData()
    A = [1 0.2 0; 0 1 0.1; 0 0 1];
    B = [1 0; 0 1; 1 -1];
    E = [1 0; 0 1; 1 1];
    Gxu = [eye(5); -eye(5); 1 -2 3 -1 2];
    Fxu = 10 + (1:size(Gxu, 1))';
    Gw = [eye(2); -eye(2)];
    Fw = [2; 3; 1; 2];
    data = {A, B, E, Gxu, Fxu, Gw, Fw};
end

function verifySingleArgumentRejections(testCase, data, index, invalid_values)
    for i = 1:size(invalid_values, 1)
        candidate = replaceArgument(data, index, invalid_values{i, 1});
        verifyRejected(testCase, candidate, invalid_values{i, 2});
    end
end

function verifyRejected(testCase, data, expected_message)
    try
        validateComputeRCISInputs(data{:});
    catch exception
        verifyTrue(testCase, contains(exception.message, expected_message), ...
            sprintf('Unexpected rejection: %s', exception.message));
        return;
    end
    verifyFail(testCase, sprintf( ...
        'Expected rejection containing message text: "%s".', expected_message));
end

function data = replaceArgument(data, index, value)
    data{index} = value;
end

function data = replaceArguments(data, indices, values)
    for i = 1:numel(indices)
        data{indices(i)} = values{i};
    end
end

function matrix = replaceEntry(matrix, index, value)
    matrix(index) = value;
end
