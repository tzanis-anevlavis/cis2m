function tests = test_construct_shrunk_safe_set_constraints
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    matlab_dir = fileparts(fileparts(mfilename('fullpath')));
    testCase.addTeardown(@path, path);
    addpath(matlab_dir, fullfile(matlab_dir, 'support_functions'));

    % The implementation uses MPT's Polyhedron and support evaluation.
    assumeTrue(testCase, exist('Polyhedron', 'class') == 8 && ...
        exist('mpt_init', 'file') == 2, 'Requires MPT3.');
    mpt_init;
end

function testSixStateJointSetUsesAnalyticBoxSupports(testCase)
    % Three Brunovsky chains and a translated, asymmetric 2D disturbance test
    % every shrinking stage without using MPT to calculate expected supports.
    chain_lengths = [3 2 1];
    [Ac, Bc] = brunovskyChains(chain_lengths);
    n = size(Ac, 1);
    m = size(Bc, 2);
    nu = max(chain_lengths);
    q = 4;
    Ec = [0 1; 1 -1; 2 0; -1 2; 0.5 1; -2 -0.5];

    mixed_Gz = [1 -2 0 1 0 3;
                0 1 -1 2 -2 0;
                2 0 1 -1 1 -2];
    Gz = [eye(n); -eye(n); mixed_Gz; zeros(m, n)];
    Gr = [zeros(2*n, m);
        1 0 -2;
        0 3 1;
        -1 2 0;
        eye(m)];
    Gc = [Gz Gr];
    Fc = 50 + (1:size(Gc, 1))';

    lower = [-1 1];
    upper = [2 3];
    Gw = [eye(2); -eye(2)];
    Fw = [upper'; -lower'];
    [G_k, F_k] = constructShrunkSafeSetConstraints( ...
        Ac, Bc, Gc, Fc, Ec, Gw, Fw, q, nu);

    % The sequence contains S_0 through S_(nu+q-1).
    N = nu + q;
    verifySize(testCase, G_k, [N 1]);
    verifySize(testCase, F_k, [N 1]);

    % Compute each box support analytically. For a box [lower, upper],
    % h_W(d) is the sum of the larger endpoint product in each coordinate.
    expected_F = Fc;
    A_curr = eye(n);
    pure_input_rows = (size(Gc, 1) - m + 1):size(Gc, 1);
    for index = 1:N
        if (index > 1 && index <= nu + 1)
            directions = Gz*A_curr*Ec;
            expected_support = sum(max(directions.*lower, directions.*upper), 2);
            expected_F = expected_F - expected_support;
            A_curr = A_curr*Ac;
        end

        % Pontryagin subtraction changes only the right-hand sides; rows with
        % zero state normal remain unchanged even if they constrain inputs.
        verifyEqual(testCase, G_k{index}, Gc);
        verifyEqual(testCase, F_k{index}, expected_F, 'AbsTol', 1e-9);
        verifyEqual(testCase, F_k{index}(pure_input_rows), Fc(pure_input_rows));
    end

    % The -e1 state normal has support h_W([0,-1]') = -1 at the first
    % stage, so a translated disturbance can legitimately loosen that row.
    verifyGreaterThan(testCase, F_k{2}(n + 1), Fc(n + 1));
    for index = (nu + 2):N
        verifyEqual(testCase, F_k{index}, F_k{nu + 1});
    end
end

function testSupportMethodMatchesMptPontryaginDifference(testCase)
    % This is the direct regression against the previous implementation.
    % MPT constructs each Pontryagin difference geometrically, while the
    % function under test updates the original H-representation by supports.
    [Ac, Bc] = brunovskyChains([2 1]);
    n = size(Ac, 1);
    m = size(Bc, 2);
    nu = 2;
    q = 3;
    Ec = [1 -0.5; -1 2; 0.75 1];

    Gc = [eye(n + m); -eye(n + m);
        1 -2 1 0.5 -1;
        -1 0 2 -2 1];
    Fc = 20 + (1:size(Gc, 1))';

    % MPT's minus implementation is a valid comparison oracle when the
    % disturbance contains the origin, as this asymmetric box does.
    lower = [-1 -0.5];
    upper = [1.5 2];
    Gw = [eye(2); -eye(2)];
    Fw = [upper'; -lower'];
    W = Polyhedron('A', Gw, 'b', Fw);

    [G_k, F_k] = constructShrunkSafeSetConstraints( ...
        Ac, Bc, Gc, Fc, Ec, Gw, Fw, q, nu);

    % Embed each state disturbance contribution in joint [z; r] space. The
    % zero input rows enforce the paper's subtraction by D x {0}, not by a
    % disturbance acting on the Brunovsky input coordinates.
    expected_set = Polyhedron('A', Gc, 'b', Fc);
    A_curr = eye(n);
    for index = 1:(nu + q)
        if (index > 1 && index <= nu + 1)
            joint_disturbance_map = [A_curr*Ec; zeros(m, size(Ec, 2))];
            expected_set = expected_set - full(joint_disturbance_map)*W;
            A_curr = A_curr*Ac;
        end

        support_set = Polyhedron('A', G_k{index}, 'b', F_k{index});
        verifyFalse(testCase, support_set.isEmptySet());
        verifyFalse(testCase, expected_set.isEmptySet());
        verifyTrue(testCase, support_set == expected_set, ...
            sprintf('Support and MPT constructions differ at cell %d.', index));
    end
end

function testNonBoxDisturbanceUsesSupportInEveryDirection(testCase)
    % W is the triangle conv{(0,0), (2,0), (0,1)}. Expected supports are
    % computed directly from these vertices, independently of MPT.
    vertices = [0 0; 2 0; 0 1];
    Gw = [-1 0; 0 -1; 0.5 1];
    Fw = [0; 0; 1];

    Ac = [0 1 0; 0 0 1; 0 0 0];
    Bc = [0; 0; 1];
    Ec = [1 0; 0 1; 1 -1];
    Gz = [1 2 -1; 1 2 -1; -2 0 3; 0 -1 2; 0 0 0];
    Gr = [2; -3; 1; 0; 4];
    Gc = [Gz Gr];
    Fc = [11; 13; 17; 19; 23];
    nu = 3;
    q = 2;

    [G_k, F_k] = constructShrunkSafeSetConstraints( ...
        Ac, Bc, Gc, Fc, Ec, Gw, Fw, q, nu);
    expected_F = boundsFromDisturbanceVertices( ...
        Ac, Gc, Fc, Ec, vertices, q, nu);

    for index = 1:(nu + q)
        verifyEqual(testCase, G_k{index}, Gc);
        verifyEqual(testCase, F_k{index}, expected_F{index}, 'AbsTol', 1e-10);

        % Rows one and two have the same state normal but different input
        % normals, so their accumulated disturbance deductions must agree.
        actual_deduction = Fc - F_k{index};
        verifyEqual(testCase, actual_deduction(1), actual_deduction(2), ...
            'AbsTol', 1e-10);
    end
end

function testLowerDimensionalDisturbanceWithSparseInputs(testCase)
    % W = conv{(1,-1), (1,2)} is a one-dimensional segment in R^2. This
    % exercises support evaluation for bounded sets with empty interior.
    vertices = [1 -1; 1 2];
    Gw = sparse([1 0; -1 0; 0 1; 0 -1]);
    Fw = sparse([1; -1; 2; 1]);

    Ac = sparse([0 1; 0 0]);
    Bc = sparse([0; 1]);
    Ec = sparse([1 -1; 2 1]);
    Gz = sparse([1 0; -1 0; 0 1; 0 -1; 0 0]);
    Gr = sparse([2; -3; 1; 4; 1]);
    Gc = [Gz Gr];
    Fc = sparse([7; 11; 13; 17; 19]);
    nu = 2;
    q = 2;

    [G_k, F_k] = constructShrunkSafeSetConstraints( ...
        Ac, Bc, Gc, Fc, Ec, Gw, Fw, q, nu);
    expected_F = boundsFromDisturbanceVertices( ...
        Ac, Gc, Fc, Ec, vertices, q, nu);

    for index = 1:(nu + q)
        verifyEqual(testCase, G_k{index}, Gc);
        verifyEqual(testCase, full(F_k{index}), expected_F{index}, ...
            'AbsTol', 1e-10);
    end

    % The fourth state normal maps through Ec to disturbance direction
    % [-2,-1], whose support on W is -1; this row is therefore loosened.
    verifyGreaterThan(testCase, F_k{2}(4), Fc(4));
end

function testNoDisturbancePreservesJointRows(testCase)
    % With no disturbance model, every returned set must be exactly S.
    Ac = [0 1; 0 0];
    Bc = [0; 1];
    G = sparse([1 2 3; -1 0 2; 0 -1 -4]);
    F = [3; 2; 1];
    [G_k, F_k] = constructShrunkSafeSetConstraints( ...
        Ac, Bc, G, F, [], [], [], 3, 2);

    verifySize(testCase, G_k, [5 1]);
    for k = 1:5
        verifyEqual(testCase, G_k{k}, G);
        verifyEqual(testCase, F_k{k}, F);
    end
end

function testEmptyJointConstraintRepresentation(testCase)
    % A 0-by-(n+m) matrix and 0-by-1 vector canonically represent full joint
    % space. The disturbance path must preserve these dimensions at all times.
    G = zeros(0, 2);
    F = zeros(0, 1);
    [G_k, F_k] = constructShrunkSafeSetConstraints( ...
        0, 1, G, F, 1, [1; -1], [2; 1], 2, 1);

    verifySize(testCase, G_k, [3 1]);
    verifySize(testCase, F_k, [3 1]);
    for index = 1:3
        verifySize(testCase, G_k{index}, [0 2]);
        verifySize(testCase, F_k{index}, [0 1]);
    end
end

function testTranslatedDisturbanceAndZeroMap(testCase)
    G = [1 0; -1 0; 0 1];
    F = [3; 3; 2];

    % For S = [-3,3] x R and W = [1,2], Pontryagin subtraction gives
    % [-4,1] x R. This is deliberately different from S + (-W) = [-5,2].
    % The support in direction -1 is -1, so the lower bound is loosened.
    [~, F_k] = constructShrunkSafeSetConstraints( ...
        0, 1, G, F, 1, [1; -1], [2; -1], 1, 1);
    verifyEqual(testCase, F_k{2}, [1; 4; 2], 'AbsTol', 1e-8);

    % A specified disturbance with Ec = 0 leaves every row unchanged.
    [G_k, F_k] = constructShrunkSafeSetConstraints( ...
        0, 1, G, F, 0, [1; -1], [2; -1], 2, 1);
    for k = 1:3
        verifyEqual(testCase, G_k{k}, G);
        verifyEqual(testCase, F_k{k}, F);
    end
end

function testShrinkingPreservesSingletonsAndInfeasibility(testCase)
    % Shrinking [-1,1] by [-1,1] collapses the state projection to {0}; the
    % implementation must retain this lower-dimensional, nonempty result.
    G = [1 0; -1 0; 0 1; 0 -1];
    F = ones(4, 1);
    [G_k, F_k] = constructShrunkSafeSetConstraints( ...
        0, 1, G, F, 1, [1; -1], [1; 1], 1, 1);
    collapsed = Polyhedron('A', G_k{2}, 'b', F_k{2});
    verifyFalse(testCase, collapsed.isEmptySet());
    verifyTrue(testCase, collapsed.contains([0; 0.5]));
    verifyFalse(testCase, collapsed.contains([0.1; 0]));
    verifyEqual(testCase, F_k{2}, [0; 0; 1; 1], 'AbsTol', 1e-8);

    % A disturbance wider than the safe state interval makes the difference
    % empty; constructing the returned inequalities must preserve that fact.
    [G_k, F_k] = constructShrunkSafeSetConstraints( ...
        0, 1, G, F, 1, [1; -1], [2; 2], 1, 1);
    infeasible = Polyhedron('A', G_k{2}, 'b', F_k{2});
    verifyTrue(testCase, infeasible.isEmptySet());
end

function testRejectsInvalidSequenceLengths(testCase)
    data = validTestData();
    invalid_values = {
        [];     % The sequence length must be scalar.
        0;      % The sequence length must be positive.
        -1;     % Negative sequence lengths are not allowed.
        1.5;    % The sequence length must be an integer.
        Inf;    % Infinite values are not allowed.
        NaN;    % NaN values are not allowed.
        1 + 1i; % Complex values are not allowed.
        [1 2];  % Vector-valued sequence lengths are not allowed.
        '2';    % Nonnumeric values are not allowed.
        true    % Logical values are not accepted as sequence lengths.
    };

    % q and nmax have the same positive-integer contract.
    verifyInvalidArgumentValues(testCase, data, 8, 'q', invalid_values);
    verifyInvalidArgumentValues(testCase, data, 9, 'nmax', invalid_values);
end

function testRejectsInvalidConstraintAttributes(testCase)
    data = validTestData();
    Gc = data{3};
    invalid_Gc = {
        repmat('x', size(Gc));     % Gc must be numeric.
        zeros([size(Gc) 2]);       % Gc must be two-dimensional.
        Gc + 1i;                   % Complex entries are not allowed.
        replaceEntry(Gc, 1, NaN);  % NaN entries are not allowed.
        replaceEntry(Gc, 1, Inf)   % Infinite entries are not allowed.
    };
    verifyInvalidArgumentValues(testCase, data, 3, 'Gc', invalid_Gc);

    Fc = data{4};
    invalid_Fc = {
        repmat('x', size(Fc));     % Fc must be numeric.
        zeros([size(Fc) 2]);       % Fc must be two-dimensional.
        Fc + 1i;                   % Complex entries are not allowed.
        replaceEntry(Fc, 1, NaN);  % NaN entries are not allowed.
        replaceEntry(Fc, 1, Inf)   % Infinite entries are not allowed.
    };
    verifyInvalidArgumentValues(testCase, data, 4, 'Fc', invalid_Fc);
end

function testRejectsInvalidSystemAndJointDimensions(testCase)
    % Ac must be square.
    verifyRejectedMessage(testCase, ...
        {zeros(2, 3), zeros(2, 1), zeros(1, 4), 1, [], [], [], 1, 1}, ...
        'Ac must be square');

    % Bc and Ac must describe the same state dimension.
    verifyRejectedMessage(testCase, ...
        {zeros(2), zeros(3, 1), zeros(1, 3), 1, [], [], [], 1, 1}, ...
        'Rows of Bc and Ac');

    % Gc must have one column for every state and input coordinate.
    verifyError(testCase, @() constructShrunkSafeSetConstraints( ...
        0, 1, ones(2, 3), [1; 1], [], [], [], 2, 1), ...
        'cis2m:constructShrunkSafeSetConstraints:InvalidJointColumns');

    % Fc must provide exactly one bound per row of Gc.
    verifyError(testCase, @() constructShrunkSafeSetConstraints( ...
        0, 1, ones(2, 2), 1, [], [], [], 2, 1), ...
        'cis2m:constructShrunkSafeSetConstraints:InvalidJointBounds');
    verifyError(testCase, @() constructShrunkSafeSetConstraints( ...
        0, 1, ones(2, 2), [1 1], [], [], [], 2, 1), ...
        'cis2m:constructShrunkSafeSetConstraints:InvalidJointBounds');

    % Ec and Ac must describe the same state dimension.
    verifyRejectedMessage(testCase, ...
        {[0 1; 0 0], [0; 1], zeros(1, 3), 1, ones(3, 1), ...
        [1; -1], [1; 1], 1, 2}, 'Rows of Ec and Ac');
end

function testRejectsNonBrunovskySystemAndWrongNmax(testCase)
    Gc = [eye(3); -eye(3)];
    Fc = ones(6, 1);

    % Bc indicates one chain of length two, but Ac omits its required shift.
    verifyError(testCase, @() constructShrunkSafeSetConstraints( ...
        zeros(2), [0; 1], Gc, Fc, [], [], [], 1, 2), ...
        'cis2m:validateBrunovskyNormalForm:InvalidStateMatrix');

    % The supplied index must equal the chain length inferred from Ac and Bc.
    verifyError(testCase, @() constructShrunkSafeSetConstraints( ...
        [0 1; 0 0], [0; 1], Gc, Fc, [], [], [], 1, 1), ...
        'cis2m:constructShrunkSafeSetConstraints:InvalidNilpotencyIndex');
end

function testRejectsInvalidDisturbanceDescriptions(testCase)
    % Exercise every inconsistent presence combination of Ec, Gw, and Fw.
    invalid_models = {
        [], [1; -1], [], 'Without disturbance';  % Gw alone is present.
        [], [], [1; 1], 'Without disturbance';   % Fw alone is present.
        [], [1; -1], [1; 1], 'Without disturbance'; % Ec alone is absent.
        1, [], [], 'With disturbance';            % Only Ec is present.
        1, [1; -1], [], 'With disturbance';       % Fw is absent.
        1, [], [1; 1], 'With disturbance'         % Gw is absent.
    };
    for index = 1:size(invalid_models, 1)
        data = {0, 1, [1 0], 1, invalid_models{index, 1}, ...
            invalid_models{index, 2}, invalid_models{index, 3}, 1, 1};
        verifyRejectedMessage(testCase, data, invalid_models{index, 4});
    end

    % Opposing inconsistent bounds describe an empty disturbance set.
    verifyError(testCase, @() constructShrunkSafeSetConstraints( ...
        0, 1, [1 0], 1, 1, [1; -1], [-1; -1], 1, 1), ...
        'cis2m:constructShrunkSafeSetConstraints:InvalidDisturbance');

    % A single upper bound leaves the scalar disturbance unbounded below.
    verifyError(testCase, @() constructShrunkSafeSetConstraints( ...
        0, 1, [1 0], 1, 1, 1, 1, 1, 1), ...
        'cis2m:constructShrunkSafeSetConstraints:InvalidDisturbance');
end

function data = validTestData()
    data = {0, 1, [1 0; -1 0], [2; 2], [], [], [], 2, 1};
end

function verifyInvalidArgumentValues( ...
        testCase, data, argument_index, argument_name, invalid_values)
    for index = 1:numel(invalid_values)
        candidate = data;
        candidate{argument_index} = invalid_values{index};
        verifyRejectedMessage(testCase, candidate, argument_name);
    end
end

function verifyRejectedMessage(testCase, data, expected_message)
    try
        constructShrunkSafeSetConstraints(data{:});
    catch exception
        verifyTrue(testCase, contains(exception.message, expected_message), ...
            sprintf('Unexpected rejection: %s', exception.message));
        return;
    end
    verifyFail(testCase, sprintf( ...
        'Expected rejection containing message text: "%s".', expected_message));
end

function matrix = replaceEntry(matrix, index, value)
    matrix(index) = value;
end

function expected_F = boundsFromDisturbanceVertices( ...
        Ac, Gc, Fc, Ec, vertices, q, nu)
    % Independently evaluate h_W(d) = max_{w in vertices(W)} d*w and sum
    % the support contributions for Ac^0 Ec W through Ac^(nu-1) Ec W.
    N = nu + q;
    n = size(Ac, 2);
    expected_F = cell(N, 1);
    expected_F{1} = full(Fc);
    accumulated_support = zeros(size(Fc, 1), 1);
    A_curr = eye(n);

    for t = 1:nu
        directions = full(Gc(:, 1:n) * A_curr * Ec);
        support = max(directions * vertices', [], 2);
        accumulated_support = accumulated_support + support;
        expected_F{t + 1} = full(Fc) - accumulated_support;
        A_curr = A_curr * Ac;
    end

    % Nilpotency makes all sets after S_nu identical to S_nu.
    for index = (nu + 2):N
        expected_F{index} = expected_F{nu + 1};
    end
end

function [Ac, Bc] = brunovskyChains(chain_lengths)
    % Construct canonical chains consecutively in state order.
    n = sum(chain_lengths);
    m = numel(chain_lengths);
    Ac = zeros(n);
    Bc = zeros(n, m);
    first = 1;
    for channel = 1:m
        last = first + chain_lengths(channel) - 1;

        % Shift each state to its successor within the current chain.
        for state = first:(last - 1)
            Ac(state, state + 1) = 1;
        end

        % The channel input acts on the final state of its chain.
        Bc(last, channel) = 1;
        first = last + 1;
    end
end
