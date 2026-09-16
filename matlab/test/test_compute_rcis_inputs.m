function tests = test_compute_rcis_inputs
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    matlab_dir = fileparts(fileparts(mfilename('fullpath')));
    testCase.addTeardown(@path, path);
    addpath(matlab_dir, fullfile(matlab_dir, 'support_functions'));

    % Scalar integrator with box constraints |x| <= 2 and |u| <= 1.
    testCase.TestData.data = {0, 1, [], [1 0; -1 0; 0 1; 0 -1], [2; 2; 1; 1], [], []};
end

function testMissingMode(testCase)
    data = testCase.TestData.data;

    % Neither computation mode is selected.
    verifyError(testCase, @() computeRCIS(data{:}), 'cis2m:computeRCIS:MissingMode');

    % Tau alone does not select single-component mode.
    verifyError(testCase, @() computeRCIS(data{:}, 'tau', 2), ...
        'cis2m:computeRCIS:MissingMode');

    % Empty lambda and hierarchy_level values both mean "not specified."
    verifyError(testCase, @() computeRCIS(data{:}, 'lambda', [], 'hierarchy_level', []), ...
        'cis2m:computeRCIS:MissingMode');
end

function testInvalidLambda(testCase)
    values = {
        0;      % Lambda must be positive.
        -1;     % Negative values are not allowed.
        1.5;    % Lambda must be an integer.
        Inf;    % Infinite values are not allowed.
        NaN;    % NaN values are not allowed.
        1 + 1i; % Complex values are not allowed.
        [1 2];  % Lambda must be scalar.
        '3';    % Nonnumeric values are not allowed.
        true    % Logical values are not accepted as numeric loop lengths.
    };
    for i = 1:numel(values)
        verifyInvalidOption(testCase, 'lambda', values{i});
    end
end

function testInvalidTau(testCase)
    values = {
        [];    % Tau must be scalar when single-component mode is selected.
        -1;    % Tau must be nonnegative.
        0.5;   % Tau must be an integer.
        Inf;   % Infinite values are not allowed.
        NaN;   % NaN values are not allowed.
        1i;    % Complex values are not allowed.
        [0 1]; % Tau must be scalar.
        '0';   % Nonnumeric values are not allowed.
        false  % Logical values are not accepted as numeric transient lengths.
    };
    for i = 1:numel(values)
        verifyInvalidOption(testCase, 'tau', values{i}, 'lambda', 1);
    end
end

function testInvalidHierarchyLevel(testCase)
    values = {
        0;     % hierarchy_level must be positive.
        -1;    % Negative values are not allowed.
        1.5;   % hierarchy_level must be an integer.
        Inf;   % Infinite values are not allowed.
        NaN;   % NaN values are not allowed.
        1i;    % Complex values are not allowed.
        [1 2]; % hierarchy_level must be scalar.
        '3';   % Nonnumeric values are not allowed.
        true   % Logical values are not accepted as numeric hierarchy levels.
    };
    for i = 1:numel(values)
        % A supplied but invalid hierarchy level must not fall back to lambda.
        verifyInvalidOption(testCase, 'hierarchy_level', values{i}, 'lambda', 1);
    end
end

function testInvalidOutputMode(testCase)
    values = {
        [];           % is_implicit must be scalar.
        -1;           % Numeric flags must be binary.
        2;            % Numeric flags must be binary.
        0.5;          % Numeric flags must be binary.
        Inf;          % Infinite values are not binary.
        NaN;          % NaN is not binary.
        1i;           % Complex values are not allowed.
        [true false]; % Logical flags must be scalar.
        'true'        % Character values are not accepted as logical flags.
    };
    for i = 1:numel(values)
        verifyInvalidOption(testCase, 'is_implicit', values{i}, 'hierarchy_level', 1);
    end
end

function testMalformedNameValueArguments(testCase)
    data = testCase.TestData.data;

    % The arguments block rejects option names outside the documented API.
    verifyError(testCase, @() computeRCIS(data{:}, 'unknown_option', 1), ...
        'MATLAB:TooManyInputs');

    % Every option name must be followed by a value.
    verifyError(testCase, @() computeRCIS(data{:}, 'lambda'), ...
        'MATLAB:TooManyInputs');
end

function testSingleComponentDefaults(testCase)
    requirePolyhedralRuntime(testCase);
    data = testCase.TestData.data;

    % The documented defaults are tau = 0 and implicit output.
    actual = computeRCIS(data{:}, 'lambda', 1);
    expected = computeRCIS(data{:}, 'lambda', 1, 'tau', 0, 'is_implicit', true);
    verifySameSets(testCase, actual, expected);

    % The lifted coordinates are [x; v], with q = tau + lambda = 1.
    verifyEqual(testCase, actual.Dim, 2);
end

function testTransientSetsSequenceLength(testCase)
    requirePolyhedralRuntime(testCase);
    data = testCase.TestData.data;
    actual = computeRCIS(data{:}, 'tau', 1, 'lambda', 2);

    % Single-component mode returns only the requested (tau, lambda) pair.
    verifyEqual(testCase, numel(actual), 1);

    % The lifted dimension is n + m*q = 1 + 1*(1 + 2) = 4.
    verifyEqual(testCase, actual.Dim, 4);
end

function testEmptyHierarchySelectsComponent(testCase)
    requirePolyhedralRuntime(testCase);
    data = testCase.TestData.data;

    % Empty hierarchy_level falls through to single-component mode.
    actual = computeRCIS(data{:}, 'lambda', 2, 'hierarchy_level', []);
    expected = computeRCIS(data{:}, 'lambda', 2);
    verifySameSets(testCase, actual, expected);
end

function testHierarchyMatchesIndividualComponents(testCase)
    requirePolyhedralRuntime(testCase);
    data = testCase.TestData.data;
    actual = computeRCIS(data{:}, 'hierarchy_level', 3);

    % Level q contains one component for each lambda = 1, ..., q.
    assertEqual(testCase, numel(actual), 3);
    for lambda = 1:3
        % Every component satisfies tau + lambda = q.
        expected = computeRCIS(data{:}, 'lambda', lambda, 'tau', 3-lambda);
        verifySameSets(testCase, actual(lambda), expected);
    end
end

function testHierarchyOverridesUnusedOptions(testCase)
    requirePolyhedralRuntime(testCase);
    data = testCase.TestData.data;
    expected = computeRCIS(data{:}, 'hierarchy_level', 2);

    % A nonempty hierarchy_level overrides otherwise valid lambda and tau.
    actual = computeRCIS(data{:}, 'lambda', 3, 'tau', 2, 'hierarchy_level', 2);
    verifySameSets(testCase, actual, expected);

    % Overridden lambda and tau values are intentionally not validated.
    actual = computeRCIS(data{:}, 'hierarchy_level', 2, 'lambda', NaN, 'tau', {'ignored'});
    verifySameSets(testCase, actual, expected);

    % Empty overridden values are also ignored in hierarchy mode.
    actual = computeRCIS(data{:}, 'lambda', [], 'tau', [], 'hierarchy_level', 2);
    verifySameSets(testCase, actual, expected);
end

function testNumericOutputFlags(testCase)
    requirePolyhedralRuntime(testCase);
    data = testCase.TestData.data;

    % Numeric 1, logical true, and the default all select implicit output.
    implicit_set = computeRCIS(data{:}, 'lambda', 1, 'is_implicit', 1);
    logical_implicit_set = computeRCIS(data{:}, 'lambda', 1, 'is_implicit', true);
    default_set = computeRCIS(data{:}, 'lambda', 1);
    verifySameSets(testCase, implicit_set, default_set);
    verifySameSets(testCase, logical_implicit_set, default_set);

    % Numeric 0 and logical false both select projection to x-space.
    explicit_set = computeRCIS(data{:}, 'lambda', 1, 'is_implicit', 0);
    logical_flag_set = computeRCIS(data{:}, 'is_implicit', false, 'lambda', 1);
    verifySameSets(testCase, explicit_set, logical_flag_set);

    % Explicit output removes the one virtual-input coordinate.
    verifyEqual(testCase, explicit_set.Dim, 1);
end

function testMixedConstraintsHaveAnalyticRCIS(testCase)
    requirePolyhedralRuntime(testCase);
    % x+ = u, |x| <= 2, |u| <= 1, x + u <= 0.5, q = 1.
    % The constant input v must also satisfy 2*v <= 0.5 at subsequent times.
    G = [1 0; -1 0; 0 1; 0 -1; 1 1];
    F = [2; 2; 1; 1; 0.5];
    [C, D] = computeRCIS(0, 1, [], G, F, [], [], 'lambda', 1);
    expected = Polyhedron('A', [G; 0 2], 'b', [F; 0.5]);
    verifyTrue(testCase, C == expected);
    verifySize(testCase, D, [2 2]);
    verifyEqual(testCase, full(D), [0 1; 0 1]);
    [Cx, Dx] = computeRCIS(0, 1, [], G, F, [], [], 'lambda', 1, 'is_implicit', false);
    expected_x = Polyhedron('A', [1; -1], 'b', [1.5; 2]);
    verifyTrue(testCase, Cx == expected_x);
    verifyEqual(testCase, Dx, D);
end

function testStateOnlyConstraintsAndOriginalFeedbackCoordinates(testCase)
    requirePolyhedralRuntime(testCase);
    % T = 1/3, Am = 2, Bm = 1. Hence u = v - (2/3)*x and x+ = 3*v.
    G = [1 0; -1 0];
    F = [3; 3];
    [C, D] = computeRCIS(2, 3, [], G, F, [], [], 'lambda', 1);
    expected = Polyhedron('A', [1 0; -1 0; 0 1; 0 -1], 'b', [3; 3; 1; 1]);
    verifyTrue(testCase, C == expected);
    verifyEqual(testCase, full(D), [0 3; 0 1], 'AbsTol', 1e-10);
end

function testMultiInputOutputCoordinatesAndJointRows(testCase)
    requirePolyhedralRuntime(testCase);
    A = [2 1; 0 3];
    B = [2 1; 1 1];
    G = [eye(4); -eye(4); 1 -2 3 -4];
    F = 10*ones(9, 1);
    [C, D] = computeRCIS(A, B, [], G, F, [], [], 'lambda', 2, 'tau', 1);
    [~, ~, ~, ~, ~, T, ~, Am, Bm] = transformToBrunovskyNormalForm(A, B, [], G, F);
    q = 3;
    H = [1 0 0 0 0 0; 0 0 0 1 0 0];
    Pbar = [0 1 0; 0 0 1; 0 1 0];
    P = blkdiag(Pbar, Pbar);
    x = [0.2; -0.1];
    v = [0.1; -0.2; 0.3; 0.4; -0.3; 0.2];
    u = Bm \ (H*v - Am*T*x);
    verifyEqual(testCase, C.Dim, 2 + 2*q);
    verifySize(testCase, D, [8 8]);
    verifyEqual(testCase, D*[x; v], [A*x + B*u; P*v], 'AbsTol', 1e-10);
    % Independently simulate the original dynamics, recovering u at each step.
    expected_A = [];
    x_map = [eye(2) zeros(2, 2*q)];
    for t = 0:3
        if (t == 0)
            sample = 1;
        else
            sample = 2 + mod(t - 1, 2);
        end
        r_map = zeros(2, 8);
        r_map(1, 2 + sample) = 1;
        r_map(2, 2 + q + sample) = 1;
        u_map = Bm \ (r_map - Am*T*x_map);
        expected_A = [expected_A; G*[x_map; u_map]];
        x_map = A*x_map + B*u_map;
    end
    expected = Polyhedron('A', expected_A, 'b', repmat(F, 4, 1));
    verifyTrue(testCase, C == expected);
end

function testRobustJointSetAndDisturbanceDynamics(testCase)
    requirePolyhedralRuntime(testCase);
    % x+ = u + w, |w| <= 0.25. For q = 1, future mixed rows require
    % 2*v + w <= 0.5, giving v <= 0.125.
    G = [1 0; -1 0; 0 1; 0 -1; 1 1];
    F = [2; 2; 1; 1; 0.5];
    [C, D] = computeRCIS(0, 1, 1, G, F, [1; -1], [0.25; 0.25], 'lambda', 1);
    expected = Polyhedron('A', [G; 0 2], 'b', [F; 0.25]);
    verifyTrue(testCase, C == expected);
    % Check every extreme disturbance through support LPs on the output set.
    pre_normal = C.A * D;
    disturbance_support = 0.25 * abs(C.A(:, 1));
    worst_next = C.support(full(pre_normal')) + disturbance_support;
    verifyLessThanOrEqual(testCase, worst_next, C.b + 1e-7);
end

function testHierarchyDynamicsAndComponents(testCase)
    requirePolyhedralRuntime(testCase);
    G = [1 0; -1 0; 0 1; 0 -1; 1 1];
    F = [2; 2; 1; 1; 0.5];
    [components, D] = computeRCIS(0, 1, [], G, F, [], [], 'hierarchy_level', 3);
    verifyEqual(testCase, numel(components), 3);
    verifyClass(testCase, D, 'cell');
    verifyEqual(testCase, numel(D), 3);
    for lambda = 1:3
        [component, Di] = computeRCIS(0, 1, [], G, F, [], [], 'lambda', lambda, 'tau', 3 - lambda);
        verifyTrue(testCase, components(lambda) == component);
        % Every hierarchy component has its own lasso generator and therefore
        % its own lifted dynamics in the corresponding cell.
        verifyEqual(testCase, D{lambda}, Di);
    end
end

function testRobustMultiInputSystemWithUnequalChains(testCase)
    requirePolyhedralRuntime(testCase);
    % This system has T = I, Bm = I, and original nilpotency index 2.
    A = [0 1 0; 1 2 0; 0 0 3];
    B = [0 0; 1 0; 0 1];
    E = [0.05; 0.1; -0.05];
    Am = [1 2 0; 0 0 3];
    closed_A = A - B*Am;
    G = [eye(5); -eye(5); 1 -2 3 4 -5];
    F = 5*ones(11, 1);
    % W = [-1, 2], tau = 1, lambda = 1, q = 2.
    [C, D] = computeRCIS(A, B, E, G, F, [1; -1], [2; 1], 'lambda', 1, 'tau', 1);
    verifyEqual(testCase, C.Dim, 7);
    verifySize(testCase, D, [7 7]);
    x_map = [eye(3) zeros(3, 4)];
    errors = zeros(3, 1);
    expected_A = [];
    expected_b = [];
    for t = 0:3
        sample = min(t + 1, 2);
        r_map = zeros(2, 7);
        r_map(1, 3 + sample) = 1;
        r_map(2, 5 + sample) = 1;
        u_map = r_map - Am*x_map;
        expected_A = [expected_A; G*[x_map; u_map]];
        % Enumerate original state/input errors under feedback, not lifted powers.
        worst_error = max(G*[errors; -Am*errors], [], 2);
        expected_b = [expected_b; F - worst_error];
        x_map = A*x_map + B*u_map;
        next_errors = closed_A*errors;
        errors = [bsxfun(@minus, next_errors, E), bsxfun(@plus, next_errors, 2*E)];
    end
    expected = Polyhedron('A', expected_A, 'b', expected_b);
    verifyTrue(testCase, C == expected);
    verifyFalse(testCase, C.isEmptySet());
    pre_normal = C.A*D;
    disturbance_normal = C.A*[E; zeros(4, 1)];
    disturbance_support = max(-disturbance_normal, 2*disturbance_normal);
    worst_next = C.support(full(pre_normal')) + disturbance_support;
    verifyLessThanOrEqual(testCase, worst_next, C.b + 1e-7);
end

function testNominalProjectionIsControlledInvariantForOriginalSystem(testCase)
    requirePolyhedralRuntime(testCase);
    [A, B, Gxu, Fxu] = doubleIntegratorData();

    % Compute the projection through the public API, then independently form
    % its one-step controlled predecessor using the original (x,u) dynamics
    % and joint constraints.
    X = computeRCIS(A, B, [], Gxu, Fxu, [], [], ...
        'tau', 1, 'lambda', 2, 'is_implicit', false);

    verifyFalse(testCase, X.isEmptySet());
    verifyProjectedSetIsRobustlyControlledInvariant( ...
        testCase, X, A, B, [], Gxu, Fxu, []);
end

function testDisturbedProjectionIsRobustlyControlledInvariantForOriginalSystem(testCase)
    requirePolyhedralRuntime(testCase);
    [A, B, Gxu, Fxu] = doubleIntegratorData();
    E = [0.05; 0.1];
    Gw = [1; -1];
    Fw = [1; 1];
    W = Polyhedron('A', Gw, 'b', Fw);

    % The predecessor uses one common control input for all disturbances and
    % subtracts the exact support of E*W from every projected-set facet.
    X = computeRCIS(A, B, E, Gxu, Fxu, Gw, Fw, ...
        'tau', 1, 'lambda', 1, 'is_implicit', false);

    verifyFalse(testCase, X.isEmptySet());
    verifyProjectedSetIsRobustlyControlledInvariant( ...
        testCase, X, A, B, E, Gxu, Fxu, W);
end

function testProjectedSetsGrowWithTransientAsInRemarkFive(testCase)
    requirePolyhedralRuntime(testCase);
    [A, B, Gxu, Fxu] = doubleIntegratorData();
    lambda = 1;
    projected_sets = cell(3, 1);

    % Remark 5(1): for fixed lambda, C_x(tau+1,lambda) contains
    % C_x(tau,lambda). Use three consecutive transient lengths.
    for tau = 0:2
        projected_sets{tau + 1} = computeRCIS( ...
            A, B, [], Gxu, Fxu, [], [], ...
            'tau', tau, 'lambda', lambda, 'is_implicit', false);
        verifyFalse(testCase, projected_sets{tau + 1}.isEmptySet());
    end

    verifyTrue(testCase, projected_sets{1} <= projected_sets{2});
    verifyTrue(testCase, projected_sets{2} <= projected_sets{3});

    % This fixture exercises strict, rather than merely equal, inclusions.
    verifyFalse(testCase, projected_sets{2} <= projected_sets{1});
    verifyFalse(testCase, projected_sets{3} <= projected_sets{2});
end

function testProjectedSetsGrowWithPeriodMultiplesAsInRemarkFive(testCase)
    requirePolyhedralRuntime(testCase);
    [A, B, Gxu, Fxu] = doubleIntegratorData();
    tau = 1;
    periods = [1 2 4];
    projected_sets = cell(size(periods));

    % Remark 5(2): for fixed tau, C_x(tau,lambda) contains
    % C_x(tau,lambda') whenever lambda is a positive integer multiple of
    % lambda'. Check the divisibility chain 1 | 2 | 4.
    for i = 1:numel(periods)
        projected_sets{i} = computeRCIS( ...
            A, B, [], Gxu, Fxu, [], [], ...
            'tau', tau, 'lambda', periods(i), 'is_implicit', false);
        verifyFalse(testCase, projected_sets{i}.isEmptySet());
    end

    verifyTrue(testCase, projected_sets{1} <= projected_sets{2});
    verifyTrue(testCase, projected_sets{2} <= projected_sets{3});

    % This fixture exercises strict, rather than merely equal, inclusions.
    verifyFalse(testCase, projected_sets{2} <= projected_sets{1});
    verifyFalse(testCase, projected_sets{3} <= projected_sets{2});
end

function testEmptyAndLowerDimensionalExplicitOutputs(testCase)
    requirePolyhedralRuntime(testCase);
    G = [1 0; -1 0; 0 1; 0 -1];
    % Fix x = u = 0; the explicit projection must retain the equality x = 0.
    C = computeRCIS(0, 1, [], G, zeros(4, 1), [], [], 'lambda', 1, 'is_implicit', false);
    verifyEqual(testCase, C.Dim, 1);
    verifyFalse(testCase, C.isEmptySet());
    verifyTrue(testCase, C.contains(0));
    verifyFalse(testCase, C.contains(1));
    % Inconsistent state constraints must remain inconsistent after lifting.
    C = computeRCIS(0, 1, [], G, [-1; -1; 1; 1], [], [], 'lambda', 1);
    verifyTrue(testCase, C.isEmptySet());
    Cx = computeRCIS(0, 1, [], G, [-1; -1; 1; 1], [], [], 'lambda', 1, 'is_implicit', false);
    verifyTrue(testCase, Cx.isEmptySet());
end

function verifyInvalidOption(testCase, name, value, varargin)
    data = testCase.TestData.data;
    try
        computeRCIS(data{:}, varargin{:}, name, value);
    catch exception
        % Distinguish option rejection from a later failure in MATLAB/MPT setup.
        verifyTrue(testCase, startsWith(exception.identifier, 'MATLAB:computeRCIS:'), ...
            exception.message);
        verifyTrue(testCase, contains(exception.message, name), exception.message);
        return;
    end
    verifyFail(testCase, ['Expected rejection of invalid option: ' name]);
end

function requirePolyhedralRuntime(testCase)
    persistent is_initialized
    assumeTrue(testCase, exist('Polyhedron', 'class') == 8 && ...
        exist('mpt_init', 'file') == 2 && exist('ctrb', 'file') == 2, ...
        'Integration tests require MPT3 and Control System Toolbox.');
    if (isempty(is_initialized))
        mpt_init;
        is_initialized = true;
    end
end

function verifySameSets(testCase, actual, expected)
    assertEqual(testCase, numel(actual), numel(expected));
    for i = 1:numel(actual)
        verifyEqual(testCase, actual(i).Dim, expected(i).Dim);
        verifyTrue(testCase, actual(i) == expected(i));
    end
end

function verifyProjectedSetIsRobustlyControlledInvariant( ...
        testCase, X, A, B, E, Gxu, Fxu, W)
    % X is an RCIS exactly when X is contained in the projection onto x of
    % all safe (x,u) pairs whose successors remain in X for every w in W.
    if (isempty(E))
        disturbance_support = zeros(size(X.b));
    else
        disturbance_directions = X.A*E;
        disturbance_support = W.support(full(disturbance_directions'));
    end

    predecessor_xu = Polyhedron( ...
        'A', [X.A*A X.A*B; Gxu], ...
        'b', [X.b - disturbance_support; Fxu]);
    predecessor = predecessor_xu.projection(1:size(A, 2), 'ifourier');

    verifyTrue(testCase, X <= predecessor, ...
        'Projected set is not robustly controlled invariant for the original system.');
end

function [A, B, Gxu, Fxu] = doubleIntegratorData()
    % Bounded joint constraints with two mixed facets produce nontrivial,
    % strictly growing projected sets for the selected Remark 5 comparisons.
    A = [1 1; 0 1];
    B = [0.5; 1];
    Gxu = [eye(3); -eye(3); 1 0.4 0.3; -0.6 1 0.5];
    Fxu = [4; 2; 1; 4; 2; 1; 2.5; 1.8];
end
