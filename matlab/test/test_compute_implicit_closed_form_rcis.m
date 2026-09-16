function tests = test_compute_implicit_closed_form_rcis
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    matlab_dir = fileparts(fileparts(mfilename('fullpath')));
    testCase.addTeardown(@path, path);
    addpath(matlab_dir, fullfile(matlab_dir, 'support_functions'));
end

function testScalarLassoImposesCurrentMixedConstraint(testCase)
    % Check every returned matrix against a hand-computed scalar example.
    % The mixed inequality z + r <= 1.5 ensures that substituting r = H * v
    % and propagating the constraint through A_lifted are both observable.
    G = [1 0; -1 0; 0 1; 0 -1; 1 1];
    F = [2; 2; 1; 1; 1.5];
    [C, f, A_lifted, H, P] = computeImplicitClosedFormRCIS(0, 1, {G; G}, {F; F}, 1, 0, 1);
    expected_next = [0 1; 0 -1; 0 1; 0 -1; 0 2];

    verifyEqual(testCase, full(C), [G; expected_next]);
    verifyEqual(testCase, f, [F; F]);
    verifyEqual(testCase, full(A_lifted), [0 1; 0 1]);
    verifyEqual(testCase, full(H), 1);
    verifyEqual(testCase, full(P), 1);
    verifyAliftedIsEventuallyPeriodic(testCase, A_lifted, 1, 0, 1);
end

function testSixStateThreeInputLiftMatchesTimeIndexedMaps(testCase)
    % Verify the complete lifted construction against independently updated
    % state and virtual-input maps. Three unequal chains, a two-step transient,
    % and a three-step period give a 21-dimensional lifted system with eight
    % constraint stages having different row counts and dense mixed normals.
    chain_lengths = [3 2 1];
    [Ac, Bc] = brunovskyChains(chain_lengths);
    n = size(Ac, 1);
    m = size(Bc, 2);
    nu = max(chain_lengths);
    tau = 2;
    lambda = 3;
    q = tau + lambda;
    N = nu + q;

    G_k = cell(N, 1);
    F_k = cell(N, 1);
    entry = 1;
    for t = 1:N
        rows = 2 + mod(t, 4);
        entries = entry:(entry + rows*(n + m) - 1);
        G_k{t} = reshape(sin(entries), rows, n + m);
        F_k{t} = 20 + t + (1:rows)'/10;
        entry = entries(end) + 1;
    end

    [C, f, A_lifted, H, P] = computeImplicitClosedFormRCIS(Ac, Bc, G_k, F_k, lambda, tau, nu);
    H_bar = [1 zeros(1, q - 1)];
    P_bar = [zeros(q - 1, 1) eye(q - 1);
            zeros(1, tau) 1 zeros(1, lambda - 1)];
    expected_H = kron(eye(m), H_bar);
    expected_P = kron(eye(m), P_bar);
    expected_A_lifted = [Ac Bc*expected_H;
                         zeros(m*q, n) expected_P];

    verifySize(testCase, C, [sum(cellfun(@(G) size(G, 1), G_k)), n + m*q]);
    verifyEqual(testCase, f, cat(1, F_k{:}));
    verifyEqual(testCase, full(H), expected_H);
    verifyEqual(testCase, full(P), expected_P);
    verifyEqual(testCase, full(A_lifted), expected_A_lifted);
    verifyEqual(testCase, full(P^(tau + lambda)), full(P^tau));
    verifyAliftedIsEventuallyPeriodic(testCase, A_lifted, nu, tau, lambda);

    % Build each constraint block directly from the time-indexed maps
    % (z_0,v_0) -> z_t and (z_0,v_0) -> v_t, without using powers of
    % A_lifted or repeating the implementation's block assembly.
    state_map = [eye(n) zeros(n, m*q)];
    v_map = [zeros(m*q, n) eye(m*q)];
    expected_blocks = cell(N, 1);
    for t = 1:N
        r_map = expected_H*v_map;
        expected_blocks{t} = G_k{t}*[state_map; r_map];
        state_map = Ac*state_map + Bc*r_map;
        v_map = expected_P*v_map;
    end
    verifyEqual(testCase, full(C), cat(1, expected_blocks{:}), 'AbsTol', 2e-13);

    % Also check the channel-major lasso samples over multiple loop traversals.
    samples = reshape(1:(m*q), q, m);
    v = samples(:);
    for t = 0:(q + 2*lambda)
        if (t < tau)
            sample = t + 1;
        else
            sample = tau + mod(t - tau, lambda) + 1;
        end
        verifyEqual(testCase, H*v, samples(sample, :)');
        v = P*v;
    end
end

function testTransientAndPeriodicEdgeCases(testCase)
    % Exercise pure-periodic, unit-period, and mixed transient-periodic
    % lassos. For each pair, simulate the lifted dynamics directly and check
    % both the finite constraint stack and periodic sample selection beyond N.
    Ac = [0 1 0; 0 0 0; 0 0 0];
    Bc = [0 0; 1 0; 0 1];
    G = [eye(5); -eye(5); 1 -2 3 4 -5];
    pairs = [0 1; 2 1; 0 3; 1 2]; % [tau, lambda]
    for pair = 1:size(pairs, 1)
        tau = pairs(pair, 1);
        lambda = pairs(pair, 2);
        q = tau + lambda;
        N = 2 + q;
        G_k = repmat({G}, N, 1);
        F_k = repmat({ones(size(G, 1), 1)}, N, 1);
        [C, ~, A_lifted, H, P] = computeImplicitClosedFormRCIS(Ac, Bc, G_k, F_k, lambda, tau, 2);
        samples = reshape(1:(2*q), q, 2);
        v = samples(:);
        z = [2; -1; 3];
        initial = [z; v];
        expected = zeros(N*size(G, 1), 1);
        for t = 0:(N + lambda - 1)
            if (t < tau)
                sample = t + 1;
            else
                sample = tau + mod(t - tau, lambda) + 1;
            end
            r = samples(sample, :)';
            verifyEqual(testCase, H*v, r);
            if (t < N)
                rows = t*size(G, 1) + (1:size(G, 1));
                expected(rows) = G*[z; r];
            end
            next_z = Ac*z + Bc*r;
            verifyEqual(testCase, A_lifted*[z; v], [next_z; P*v]);
            z = next_z;
            v = P*v;
        end
        verifyEqual(testCase, C*initial, expected);
        verifyEqual(testCase, full(P^(tau + lambda)), full(P^tau));
        verifyAliftedIsEventuallyPeriodic(testCase, A_lifted, 2, tau, lambda);
    end
end

function testLiftedSetIsPositivelyInvariantForLassoEdgeCases(testCase)
    requireMpt(testCase);

    % In the disturbance-free case, verify the defining PIS property from
    % Theorem 3.5: C_lifted is contained in its one-step predecessor under
    % A_lifted. Use three unequal Brunovsky chains, mixed state-input
    % constraints, and four lasso shapes covering zero transient, unit period,
    % and nontrivial transient-period combinations.
    chain_lengths = [3 2 1];
    [Ac, Bc] = brunovskyChains(chain_lengths);
    n = size(Ac, 1);
    m = size(Bc, 2);
    nu = max(chain_lengths);
    safe_bounds = [4; 3.5; 3; 2.5; 2; 1.5; 1.2; 1; 0.8];
    mixed_rows = [1 -0.5 0.25 0 0.5 -0.25 0.75 -0.5 0.25;
                 -0.25 0.75 0 -0.5 0.25 1 -0.5 0.25 0.5];
    Gc = [eye(n + m); -eye(n + m); mixed_rows];
    Fc = [safe_bounds; safe_bounds; 0.7*abs(mixed_rows)*safe_bounds];
    lasso_shapes = [0 1; 0 3; 2 1; 2 3]; % [tau, lambda]

    for index = 1:size(lasso_shapes, 1)
        tau = lasso_shapes(index, 1);
        lambda = lasso_shapes(index, 2);
        q = tau + lambda;
        [G_k, F_k] = constructShrunkSafeSetConstraints(Ac, Bc, Gc, Fc, [], [], [], q, nu);
        [C, f, A_lifted] = computeImplicitClosedFormRCIS(Ac, Bc, G_k, F_k, lambda, tau, nu);
        lifted_set = Polyhedron('A', C, 'b', f);
        diagnostic = sprintf( ...
            'Lifted set is not positively invariant for tau=%d, lambda=%d.', ...
            tau, lambda);

        verifyFalse(testCase, lifted_set.isEmptySet(), diagnostic);
        verifyTrue(testCase, lifted_set.isBounded(), diagnostic);
        verifyTrue(testCase, isPositivelyInvariant(lifted_set, A_lifted), diagnostic);
        verifyAliftedIsEventuallyPeriodic(testCase, A_lifted, nu, tau, lambda);
    end
end

function testLiftedSetIsRobustlyPositivelyInvariant(testCase)
    requireMpt(testCase);

    % In the disturbed case, verify the defining RPIS property from Theorem
    % 3.5 by containment in the one-step robust predecessor. Exercise two
    % Brunovsky chains, mixed state-input constraints, a lasso with a transient
    % and a nontrivial period, and two non-box disturbance sets: one containing
    % zero and one translated away from zero.
    chain_lengths = [2 1];
    [Ac, Bc] = brunovskyChains(chain_lengths);
    n = size(Ac, 1);
    m = size(Bc, 2);
    nu = max(chain_lengths);
    tau = 1;
    lambda = 2;
    q = tau + lambda;

    safe_bounds = [4; 3; 2; 1.5; 1.25];
    mixed_rows = [1 -0.5 0.25 0.75 -0.5;
                 -0.25 1 -0.5 -0.25 0.75];
    Gc = [eye(n + m); -eye(n + m); mixed_rows];
    Fc = [safe_bounds; safe_bounds; 0.75*abs(mixed_rows)*safe_bounds];
    Ec = [1 0.25; -0.5 0.75; 0.25 1];
    Gw = [eye(2); -eye(2); 1 1];
    disturbance_bounds = {
        [0.10; 0.08; 0.12; 0.09; 0.12];  % Non-box W containing the origin.
        [0.12; 0.08; -0.04; 0.02; 0.16]  % Translated non-box W excluding zero.
    };

    for disturbance = 1:numel(disturbance_bounds)
        Fw = disturbance_bounds{disturbance};
        W = Polyhedron('A', Gw, 'b', Fw);
        [G_k, F_k] = constructShrunkSafeSetConstraints(Ac, Bc, Gc, Fc, Ec, Gw, Fw, q, nu);
        [C, f, A_lifted] = computeImplicitClosedFormRCIS(Ac, Bc, G_k, F_k, lambda, tau, nu);
        lifted_set = Polyhedron('A', C, 'b', f);
        E_lifted = [Ec; zeros(m*q, size(Ec, 2))];

        diagnostic = sprintf( ...
            'Lifted set failed for disturbance case %d.', disturbance);
        verifyFalse(testCase, lifted_set.isEmptySet(), diagnostic);
        verifyTrue(testCase, lifted_set.isBounded(), diagnostic);
        verifyRobustPositiveInvariance(testCase, lifted_set, A_lifted, E_lifted, W, diagnostic);
        verifyAliftedIsEventuallyPeriodic(testCase, A_lifted, nu, tau, lambda);
    end
end

function testRandomNominalLiftedSetsArePositivelyInvariant(testCase)
    requireMpt(testCase);

    % Exercise ten reproducible nominal systems with n = q = 1, ..., 10.
    % Random chain partitions, mixed safe-set facets, and lasso shapes cover
    % interactions not represented by the fixed examples. Explicit
    % nonemptiness prevents a vacuous positive-invariance result.
    seed = 20260916;
    previous_rng = rng;
    testCase.addTeardown(@() rng(previous_rng));
    rng(seed, 'twister');

    for trial = 1:10
        n = trial;
        q = trial;
        m = min(n, 1 + mod(trial - 1, 4));
        chain_lengths = randomChainLengths(n, m);
        [Ac, Bc] = brunovskyChains(chain_lengths);
        nu = max(chain_lengths);
        if (mod(trial, 2) == 1)
            lambda = q;
        else
            lambda = q/2;
        end
        tau = q - lambda;
        [Gc, Fc, state_bounds, input_bounds] = randomJointSafeSet(n, m);

        [G_k, F_k] = constructShrunkSafeSetConstraints(Ac, Bc, Gc, Fc, [], [], [], q, nu);
        [C, f, A_lifted] = computeImplicitClosedFormRCIS(Ac, Bc, G_k, F_k, lambda, tau, nu);
        lifted_set = Polyhedron('A', C, 'b', f);
        diagnostic = randomTrialDiagnostic('nominal', seed, trial, n, m, tau, lambda);

        verifyFalse(testCase, lifted_set.isEmptySet(), diagnostic);
        verifyLiftedSetIsContainedInBoundingBox( ...
            testCase, lifted_set, state_bounds, input_bounds, q, diagnostic);
        verifyTrue(testCase, isPositivelyInvariant(lifted_set, A_lifted), diagnostic);
        verifyAliftedIsEventuallyPeriodic(testCase, A_lifted, nu, tau, lambda);
    end
end

function testRandomDisturbedLiftedSetsAreRobustlyPositivelyInvariant(testCase)
    requireMpt(testCase);

    % Exercise ten reproducible disturbed systems with n = q = 1, ..., 10.
    % Each trial uses a small random disturbance map and a bounded non-box W.
    % The robust predecessor is computed independently from the returned set.
    seed = 20260917;
    previous_rng = rng;
    testCase.addTeardown(@() rng(previous_rng));
    rng(seed, 'twister');

    for trial = 1:10
        n = trial;
        q = trial;
        m = min(n, 1 + mod(trial + 1, 4));
        chain_lengths = randomChainLengths(n, m);
        [Ac, Bc] = brunovskyChains(chain_lengths);
        nu = max(chain_lengths);
        if (mod(trial, 2) == 1)
            lambda = q;
        else
            lambda = q/2;
        end
        tau = q - lambda;
        [Gc, Fc, state_bounds, input_bounds] = randomJointSafeSet(n, m);

        disturbance_dimension = min(3, n);
        Ec = 0.05 * randn(n, disturbance_dimension) / sqrt(disturbance_dimension);
        [Gw, Fw] = randomDisturbanceSet(disturbance_dimension);
        W = Polyhedron('A', Gw, 'b', Fw);

        [G_k, F_k] = constructShrunkSafeSetConstraints(Ac, Bc, Gc, Fc, Ec, Gw, Fw, q, nu);
        [C, f, A_lifted] = computeImplicitClosedFormRCIS(Ac, Bc, G_k, F_k, lambda, tau, nu);
        lifted_set = Polyhedron('A', C, 'b', f);
        E_lifted = [Ec; zeros(m*q, disturbance_dimension)];
        diagnostic = randomTrialDiagnostic('disturbed', seed, trial, n, m, tau, lambda);

        verifyFalse(testCase, lifted_set.isEmptySet(), diagnostic);
        verifyLiftedSetIsContainedInBoundingBox( ...
            testCase, lifted_set, state_bounds, input_bounds, q, diagnostic);
        verifyRobustPositiveInvariance(testCase, lifted_set, A_lifted, E_lifted, W, diagnostic);
        verifyAliftedIsEventuallyPeriodic(testCase, A_lifted, nu, tau, lambda);
    end
end

function testVariableRowsEmptyBlocksAndInvalidWidth(testCase)
    % Constraint stages may have different row counts, including no rows.
    % Check mixed empty/nonempty stages, the canonical full-space output, and
    % rejection of a stage having the wrong ambient dimension.
    G_k = {zeros(0, 2); [1 2]; [0 1; 0 -1]};
    F_k = {zeros(0, 1); 3; [1; 1]};
    [C, f, A_lifted] = computeImplicitClosedFormRCIS(0, 1, G_k, F_k, 2, 0, 1);
    verifyEqual(testCase, full(C), [0 1 2; 0 1 0; 0 -1 0]);
    verifyEqual(testCase, f, [3; 1; 1]);
    verifyAliftedIsEventuallyPeriodic(testCase, A_lifted, 1, 0, 2);

    [C, f, A_lifted] = computeImplicitClosedFormRCIS(0, 1, ...
        repmat({zeros(0, 2)}, 3, 1), repmat({zeros(0, 1)}, 3, 1), 2, 0, 1);
    verifySize(testCase, C, [0 3]);
    verifySize(testCase, f, [0 1]);
    verifyAliftedIsEventuallyPeriodic(testCase, A_lifted, 1, 0, 2);

    verifyError(testCase, @() computeImplicitClosedFormRCIS(0, 1, ...
        repmat({ones(1, 3)}, 3, 1), {1; 1; 1}, 2, 0, 1), ...
        'cis2m:computeImplicitClosedFormRCIS:InvalidConstraintMatrix');
end

function testAcceptsSparseMatricesAndRowCellArrays(testCase)
    % Verify that sparse system/constraint matrices and row-oriented cell
    % arrays are accepted without changing dimensions or numerical values.
    G = sparse([1 0; -1 0; 0 1; 0 -1]);
    F = sparse(ones(4, 1));

    % Cell orientation is immaterial because stages are indexed linearly.
    [C, f, A_lifted, H, P] = computeImplicitClosedFormRCIS(sparse(0), sparse(1), {G, G}, {F, F}, 1, 0, 1);

    verifySize(testCase, C, [8 2]);
    verifyEqual(testCase, f, [F; F]);
    verifyEqual(testCase, full(A_lifted), [0 1; 0 1]);
    verifyEqual(testCase, full(H), 1);
    verifyEqual(testCase, full(P), 1);
    verifyAliftedIsEventuallyPeriodic(testCase, A_lifted, 1, 0, 1);
end

function testRejectsInvalidSystemMatrices(testCase)
    % Change exactly one system-matrix property at a time and verify that
    % numeric type, array dimensionality, reality, finiteness, and shape.
    data = validTestData();
    A = data{1};
    B = data{2};

    invalid_A = {
        'x', 'Ac';                  % Ac must be numeric.
        zeros(1, 1, 2), 'Ac';       % Ac must be two-dimensional.
        A + 1i, 'Ac';               % Ac must be real.
        NaN, 'Ac';                  % Ac cannot contain NaN.
        Inf, 'Ac';                  % Ac cannot contain infinity.
        zeros(1, 2), 'square'       % Ac must be square.
    };
    verifyInvalidArgumentValues(testCase, data, 1, invalid_A);

    invalid_B = {
        'x', 'Bc';                  % Bc must be numeric.
        zeros(1, 1, 2), 'Bc';       % Bc must be two-dimensional.
        B + 1i, 'Bc';               % Bc must be real.
        NaN, 'Bc';                  % Bc cannot contain NaN.
        Inf, 'Bc';                  % Bc cannot contain infinity.
        zeros(2, 1), 'same number'  % Bc and Ac must have equal row counts.
    };
    verifyInvalidArgumentValues(testCase, data, 2, invalid_B);
end

function testRejectsInvalidLassoParameters(testCase)
    % Change exactly one scalar-parameter property at a time. Together these
    % cases exercise the distinct domains lambda > 0, tau >= 0, and nmax > 0.
    data = validTestData();

    invalid_lambda = {
        'x', 'lambda';       % lambda must be numeric.
        [1 1], 'lambda';     % lambda must be scalar.
        1 + 1i, 'lambda';    % lambda must be real.
        NaN, 'lambda';       % lambda cannot be NaN.
        Inf, 'lambda';       % lambda cannot be infinite.
        0, 'lambda';         % lambda must be positive.
        -1, 'lambda';        % lambda cannot be negative.
        1.5, 'lambda'        % lambda must be an integer.
    };
    verifyInvalidArgumentValues(testCase, data, 5, invalid_lambda);

    invalid_tau = {
        'x', 'tau';          % tau must be numeric.
        [0 1], 'tau';        % tau must be scalar.
        1i, 'tau';           % tau must be real.
        NaN, 'tau';          % tau cannot be NaN.
        Inf, 'tau';          % tau cannot be infinite.
        -1, 'tau';           % tau must be nonnegative.
        0.5, 'tau'           % tau must be an integer.
    };
    verifyInvalidArgumentValues(testCase, data, 6, invalid_tau);

    invalid_nmax = {
        'x', 'nmax';         % nmax must be numeric.
        [1 1], 'nmax';       % nmax must be scalar.
        1 + 1i, 'nmax';      % nmax must be real.
        NaN, 'nmax';         % nmax cannot be NaN.
        Inf, 'nmax';         % nmax cannot be infinite.
        0, 'nmax';           % nmax must be positive.
        -1, 'nmax';          % nmax cannot be negative.
        1.5, 'nmax'          % nmax must be an integer.
    };
    verifyInvalidArgumentValues(testCase, data, 7, invalid_nmax);
end

function testRejectsInvalidConstraintEntries(testCase)
    % Corrupt one stage at a time and verify every matrix/vector predicate:
    % numeric type, array dimensionality, reality, finiteness, orientation,
    % and dimensions.
    data = validTestData();
    G = data{3}{1};
    F = data{4}{1};

    invalid_G = {
        repmat('x', size(G)), 'G_k{1}';       % Constraint matrix must be numeric.
        zeros([size(G) 2]), 'G_k{1}';         % Constraint matrix must be 2D.
        G + 1i, 'G_k{1}';                     % Constraint matrix must be real.
        replaceEntry(G, 1, NaN), 'G_k{1}';    % Constraint matrix cannot contain NaN.
        replaceEntry(G, 1, Inf), 'G_k{1}';    % Constraint matrix cannot contain infinity.
        ones(size(G, 1), 3), 'G_k{1}'         % Matrix must have n + m columns.
    };
    for index = 1:size(invalid_G, 1)
        candidate = data;
        candidate{3}{1} = invalid_G{index, 1};
        verifyRejectedMessage(testCase, candidate, invalid_G{index, 2});
    end

    invalid_F = {
        repmat('x', size(F)), 'F_k{1}';       % Constraint vector must be numeric.
        zeros([size(F) 2]), 'F_k{1}';         % Constraint vector must be 2D.
        F + 1i, 'F_k{1}';                     % Constraint vector must be real.
        replaceEntry(F, 1, NaN), 'F_k{1}';    % Constraint vector cannot contain NaN.
        replaceEntry(F, 1, Inf), 'F_k{1}';    % Constraint vector cannot contain infinity.
        F', 'F_k{1}';                         % Constraint vector must be a column.
        F(1:(end - 1)), 'F_k{1}'              % One bound is required per matrix row.
    };
    for index = 1:size(invalid_F, 1)
        candidate = data;
        candidate{4}{1} = invalid_F{index, 1};
        verifyRejectedMessage(testCase, candidate, invalid_F{index, 2});
    end
end

function testInvalidCollectionSizesAndBounds(testCase)
    % Collection-level validation is separate from validation of each cell:
    % both arguments must be cells, both need N stages, and corresponding
    % matrices and bounds must have compatible row counts.
    G = [1 0];

    % Reject each non-cell collection independently.
    verifyError(testCase, @() computeImplicitClosedFormRCIS( ...
        0, 1, G, {1; 1}, 1, 0, 1), ...
        'cis2m:computeImplicitClosedFormRCIS:InvalidConstraintCollections');
    verifyError(testCase, @() computeImplicitClosedFormRCIS( ...
        0, 1, {G; G}, [1; 1], 1, 0, 1), ...
        'cis2m:computeImplicitClosedFormRCIS:InvalidConstraintCollections');

    % Reject each incorrect collection length independently.
    verifyError(testCase, @() computeImplicitClosedFormRCIS( ...
        0, 1, {G}, {1; 1}, 1, 0, 1), ...
        'cis2m:computeImplicitClosedFormRCIS:InvalidConstraintCount');
    verifyError(testCase, @() computeImplicitClosedFormRCIS( ...
        0, 1, {G; G}, {1}, 1, 0, 1), ...
        'cis2m:computeImplicitClosedFormRCIS:InvalidConstraintCount');

    % Reject a row vector where a column of bounds is required.
    verifyError(testCase, @() computeImplicitClosedFormRCIS( ...
        0, 1, {G; G}, {1; [1 2]}, 1, 0, 1), ...
        'cis2m:computeImplicitClosedFormRCIS:InvalidConstraintVector');
end

function testRejectsNonBrunovskySystemAndWrongNmax(testCase)
    % The closed form relies on exact Brunovsky chain dynamics and on nmax
    % matching the largest chain length; reject violations independently.
    G = [eye(3); -eye(3)];
    G_k = repmat({G}, 3, 1);
    F_k = repmat({ones(6, 1)}, 3, 1);

    % Bc indicates one chain of length two, but Ac omits its required shift.
    verifyError(testCase, @() computeImplicitClosedFormRCIS( ...
        zeros(2), [0; 1], G_k, F_k, 1, 0, 2), ...
        'cis2m:validateBrunovskyNormalForm:InvalidStateMatrix');

    % The supplied index must equal the chain length inferred from Ac and Bc.
    verifyError(testCase, @() computeImplicitClosedFormRCIS( ...
        [0 1; 0 0], [0; 1], {G; G}, {ones(6, 1); ones(6, 1)}, 1, 0, 1), ...
        'cis2m:computeImplicitClosedFormRCIS:InvalidNilpotencyIndex');
end

function data = validTestData()
    G = [1 0; -1 0; 0 1; 0 -1];
    F = ones(4, 1);
    data = {0, 1, {G; G}, {F; F}, 1, 0, 1};
end

function verifyInvalidArgumentValues(testCase, data, argument_index, invalid_values)
    for index = 1:size(invalid_values, 1)
        candidate = data;
        candidate{argument_index} = invalid_values{index, 1};
        verifyRejectedMessage(testCase, candidate, invalid_values{index, 2});
    end
end

function verifyRejectedMessage(testCase, data, expected_message)
    try
        computeImplicitClosedFormRCIS(data{:});
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

function chain_lengths = randomChainLengths(n, m)
    % Return a random composition of n into m positive chain lengths.
    if (m == 1)
        chain_lengths = n;
        return;
    end
    cuts = sort(randperm(n - 1, m - 1));
    chain_lengths = diff([0 cuts n]);
end

function [Gc, Fc, state_bounds, input_bounds] = randomJointSafeSet(n, m)
    % Box facets guarantee boundedness; normalized mixed facets cut random
    % corners while retaining the origin strictly in the safe-set interior.
    dimension = n + m;
    state_bounds = 4 + rand(n, 1);
    input_bounds = 2 + rand(m, 1);
    coordinate_bounds = [state_bounds; input_bounds];
    mixed_rows = normalizedRandomRows(min(6, dimension), dimension);
    mixed_bounds = (0.65 + 0.25*rand(size(mixed_rows, 1), 1)) ...
        .* (abs(mixed_rows)*coordinate_bounds);
    Gc = [eye(dimension); -eye(dimension); mixed_rows];
    Fc = [coordinate_bounds; coordinate_bounds; mixed_bounds];
end

function [Gw, Fw] = randomDisturbanceSet(dimension)
    % Construct a full-dimensional non-box polytope containing the origin.
    coordinate_bounds = 0.25 + 0.25*rand(dimension, 1);
    mixed_rows = normalizedRandomRows(dimension + 1, dimension);
    mixed_bounds = (0.65 + 0.25*rand(size(mixed_rows, 1), 1)) ...
        .* (abs(mixed_rows)*coordinate_bounds);
    Gw = [eye(dimension); -eye(dimension); mixed_rows];
    Fw = [coordinate_bounds; coordinate_bounds; mixed_bounds];
end

function rows = normalizedRandomRows(row_count, column_count)
    rows = randn(row_count, column_count);
    rows = rows ./ vecnorm(rows, 2, 2);
end

function diagnostic = randomTrialDiagnostic( ...
        kind, seed, trial, n, m, tau, lambda)
    diagnostic = sprintf( ...
        ['Random %s trial %d failed (seed=%d, n=%d, m=%d, ' ...
         'tau=%d, lambda=%d).'], ...
        kind, trial, seed, n, m, tau, lambda);
end

function requireMpt(testCase)
    assumeTrue(testCase, exist('Polyhedron', 'class') == 8 && ...
        exist('mpt_init', 'file') == 2, 'Requires MPT3.');
    mpt_init;
end

function verifyRobustPositiveInvariance(testCase, X, A, E, W, diagnostic)
    % X is robustly positively invariant exactly when it is contained in its
    % robust predecessor {xi | X.A*A*xi <= X.b - h_W((X.A*E)')}.
    directions = X.A*E;
    support = W.support(full(directions'));
    robust_predecessor = Polyhedron('A', X.A*A, 'b', X.b - support);
    verifyTrue(testCase, X <= robust_predecessor, diagnostic);
end

function verifyLiftedSetIsContainedInBoundingBox( ...
        testCase, lifted_set, state_bounds, input_bounds, q, diagnostic)
    % The t = 0 state constraints bound z. For each input channel ell, the
    % constraints over t = 0, ..., q - 1 bound every lasso coordinate as
    % -input_bounds(ell) <= v(ell,j) <= input_bounds(ell), j = 1, ..., q.
    virtual_input_bounds = kron(input_bounds, ones(q, 1));
    lifted_bounds = [state_bounds; virtual_input_bounds];
    dimension = numel(lifted_bounds);
    bounding_box = Polyhedron( ...
        'A', [speye(dimension); -speye(dimension)], ...
        'b', [lifted_bounds; lifted_bounds]);

    boundedness_diagnostic = sprintf( ...
        '%s Expected lifted bounding box is not bounded.', diagnostic);
    verifyTrue(testCase, bounding_box.isBounded(), boundedness_diagnostic);

    containment_diagnostic = sprintf( ...
        '%s Lifted set is not contained in its expected bounding box.', ...
        diagnostic);
    verifyTrue(testCase, lifted_set <= bounding_box, containment_diagnostic);
end

function verifyAliftedIsEventuallyPeriodic(testCase, A_lifted, nmax, tau, lambda)
    % Nilpotency of Ac removes the initial-state response after nmax steps,
    % while the lasso generator repeats after tau steps with period lambda.
    % Hence A_lifted^(t + lambda) = A_lifted^t for t >= nmax + tau.
    periodicity_start = nmax + tau;
    A_at_start = A_lifted^periodicity_start;
    A_one_period_later = A_lifted^(periodicity_start + lambda);

    verifyEqual(testCase, full(A_one_period_later), full(A_at_start), ...
        'AbsTol', 1e-13);

    % The nonzero repeating power also proves that A_lifted is not nilpotent.
    verifyGreaterThan(testCase, norm(A_at_start, 'fro'), 0);
end
