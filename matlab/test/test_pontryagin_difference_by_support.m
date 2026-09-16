function tests = test_pontryagin_difference_by_support
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    matlab_dir = fileparts(fileparts(mfilename('fullpath')));
    testCase.addTeardown(@path, path);
    addpath(fullfile(matlab_dir, 'support_functions'));
    assumeTrue(testCase, exist('Polyhedron', 'class') == 8 && ...
        exist('mpt_init', 'file') == 2, 'Requires MPT3.');
    mpt_init;
end

function testKnownTranslatedIntervalDifference(testCase)
    % For S = [-3,3] and W = [1,2], S - W = [-4,1]. This explicitly
    % distinguishes Pontryagin subtraction from S + (-W) = [-5,2].
    G = [1; -1];
    F = [3; 3];
    W = Polyhedron('A', [1; -1], 'b', [2; -1]);

    [G_difference, F_difference] = pontryaginDifferenceBySupport(G, F, 1, W);

    verifyEqual(testCase, G_difference, G);
    verifyEqual(testCase, F_difference, [1; 4], 'AbsTol', 1e-12);

    actual = Polyhedron('A', G_difference, 'b', F_difference);
    expected = Polyhedron('lb', -4, 'ub', 1);
    verifyTrue(testCase, actual == expected);

    % Do not use MPT minus as the oracle here. Its implementation retains
    % the original constraints of S, which is only valid when 0 belongs to W.
end

function testMappedTriangleMatchesVerticesAndMpt(testCase)
    % Use a non-box W and a rank-deficient map into a four-dimensional set.
    % Expected supports come from the known triangle vertices, not from MPT.
    vertices = [0 0; 2 0; 0 1];
    W = Polyhedron('A', [-1 0; 0 -1; 0.5 1], 'b', [0; 0; 1]);
    disturbance_map = [1 -1; 2 0.5; -1 2; 0 0];
    G = [eye(4); -eye(4);
        1 -2 3 4;
        -1 0.5 2 -3];
    F = 20 + (1:size(G, 1))';

    [G_difference, F_difference] = pontryaginDifferenceBySupport(G, F, disturbance_map, W);
    directions = G*disturbance_map;
    expected_support = max(directions*vertices', [], 2);

    verifyEqual(testCase, G_difference, G);
    verifyEqual(testCase, F_difference, F - expected_support, 'AbsTol', 1e-12);

    actual = Polyhedron('A', G_difference, 'b', F_difference);
    expected = Polyhedron('A', G, 'b', F) - full(disturbance_map)*W;
    verifyTrue(testCase, actual == expected);

    % The final coordinate is outside the image of W, so its pure rows have
    % zero support and retain their original bounds.
    verifyEqual(testCase, F_difference([4 8]), F([4 8]));
end

function testOneHundredRandomDifferencesMatchMpt(testCase)
    % MPT minus retains the original set constraints, so it is a valid oracle
    % for these trials because every generated W contains the origin.
    seed = 20260915;
    previous_rng = rng;
    testCase.addTeardown(@() rng(previous_rng));
    rng(seed, 'twister');

    for trial = 1:100
        state_dimension = 2 + mod(trial - 1, 9);
        if (state_dimension == 10)
            % Exercise genuinely high-dimensional P and W without requiring
            % a projection of W into a different-dimensional ambient space.
            disturbance_dimension = 10;
        else
            disturbance_dimension = 1 + mod(floor((trial - 1)/9), min(4, state_dimension));
        end

        % Box facets guarantee a bounded full-dimensional safe set. Random
        % normalized facets make its H-representation non-axis-aligned.
        state_lower = -(2 + rand(state_dimension, 1));
        state_upper = 2 + rand(state_dimension, 1);
        [extra_G, extra_F] = randomBoxCornerCuts( ...
            state_lower, state_upper, state_dimension + 2);
        G = [eye(state_dimension); -eye(state_dimension); extra_G];
        F = [state_upper; -state_lower; extra_F];
        safe_set = Polyhedron('A', G, 'b', F);

        % Positive right-hand sides place the origin strictly inside W, while
        % box facets guarantee boundedness independently of random facets.
        disturbance_lower = -(0.25 + 0.5*rand(disturbance_dimension, 1));
        disturbance_upper = 0.25 + 0.5*rand(disturbance_dimension, 1);
        [extra_Gw, extra_Fw] = randomBoxCornerCuts( ...
            disturbance_lower, disturbance_upper, ...
            disturbance_dimension + 2);
        Gw = [eye(disturbance_dimension); ...
            -eye(disturbance_dimension); extra_Gw];
        Fw = [disturbance_upper; -disturbance_lower; extra_Fw];
        W = Polyhedron('A', Gw, 'b', Fw);

        % A small rectangular map gives a nonempty difference and exercises
        % full-rank and rank-deficient images across varying dimensions.
        if (state_dimension == 10)
            disturbance_map = diag(0.1 + 0.1*rand(state_dimension, 1));
        else
            disturbance_map = 0.2*randn( ...
                state_dimension, disturbance_dimension) ...
                /sqrt(disturbance_dimension);
        end
        if (state_dimension < 10 && mod(trial, 5) == 0)
            disturbance_map(end, :) = 0;
        end

        [G_difference, F_difference] = pontryaginDifferenceBySupport( ...
            G, F, disturbance_map, W);
        support_difference = Polyhedron( ...
            'A', G_difference, 'b', F_difference);
        mpt_difference = safe_set - full(disturbance_map)*W;

        diagnostic = sprintf('Random trial %d failed with seed %d.', ...
            trial, seed);
        verifyFalse(testCase, support_difference.isEmptySet(), diagnostic);
        verifyFalse(testCase, mpt_difference.isEmptySet(), diagnostic);
        verifyTrue(testCase, support_difference == mpt_difference, diagnostic);
    end
end

function testLowerDimensionalDisturbanceMatchesVertexSupports(testCase)
    % W = conv{(1,-1), (1,2)} is bounded but has empty interior in R^2.
    vertices = [1 -1; 1 2];
    W = Polyhedron('A', [1 0; -1 0; 0 1; 0 -1], ...
        'b', [1; -1; 2; 1]);
    disturbance_map = sparse([1 -1; 2 1; -1 -2]);
    G = sparse([eye(3); -eye(3); 1 -2 3]);
    F = sparse([7; 11; 13; 17; 19; 23; 29]);

    [G_difference, F_difference] = pontryaginDifferenceBySupport(G, F, disturbance_map, W);
    expected_support = max(full(G*disturbance_map)*vertices', [], 2);

    verifyEqual(testCase, G_difference, G);
    verifyEqual(testCase, full(F_difference), full(F) - expected_support, 'AbsTol', 1e-12);

    % This translated segment does not contain the origin, so explicit vertex
    % supports are the independent oracle rather than MPT's minus overload.
end

function testZeroDirectionsAndEmptyConstraintRepresentation(testCase)
    W = Polyhedron('lb', -2, 'ub', 3);
    G = [1 0 0; -1 0 0; 0 1 0; 0 0 -1];
    F = [5; 7; 11; 13];
    disturbance_map = [2; 0; 0];

    [G_difference, F_difference] = pontryaginDifferenceBySupport(G, F, disturbance_map, W);
    verifyEqual(testCase, G_difference, G);
    verifyEqual(testCase, F_difference, [-1; 3; 11; 13], 'AbsTol', 1e-12);

    % A zero map leaves every inequality unchanged.
    [G_difference, F_difference] = pontryaginDifferenceBySupport(G, F, zeros(3, 1), W);
    verifyEqual(testCase, G_difference, G);
    verifyEqual(testCase, F_difference, F);

    % Canonical full-space H-representations retain their dimensions.
    [G_difference, F_difference] = pontryaginDifferenceBySupport(zeros(0, 3), zeros(0, 1), disturbance_map, W);
    verifySize(testCase, G_difference, [0 3]);
    verifySize(testCase, F_difference, [0 1]);
end

function testDifferenceCanBeLowerDimensionalOrEmpty(testCase)
    G = [1; -1];
    F = [1; 1];

    % Subtracting [-1,1] from itself produces the singleton {0}.
    W = Polyhedron('lb', -1, 'ub', 1);
    [G_singleton, F_singleton] = pontryaginDifferenceBySupport(G, F, 1, W);
    singleton = Polyhedron('A', G_singleton, 'b', F_singleton);
    verifyFalse(testCase, singleton.isEmptySet());
    verifyTrue(testCase, singleton.contains(0));
    verifyFalse(testCase, singleton.contains(0.1));

    % A wider disturbance produces an empty Pontryagin difference.
    W = Polyhedron('lb', -2, 'ub', 2);
    [G_empty, F_empty] = pontryaginDifferenceBySupport(G, F, 1, W);
    verifyTrue(testCase, Polyhedron('A', G_empty, 'b', F_empty).isEmptySet());
end

function testRejectsInvalidRepresentationsAndMaps(testCase)
    W = Polyhedron('lb', -1, 'ub', 1);

    % F must be a column vector with one entry per inequality.
    verifyError(testCase, @() pontryaginDifferenceBySupport( ...
        ones(2, 2), 1, ones(2, 1), W), ...
        'cis2m:pontryaginDifferenceBySupport:InvalidBounds');
    verifyError(testCase, @() pontryaginDifferenceBySupport( ...
        ones(2, 2), [1 1], ones(2, 1), W), ...
        'cis2m:pontryaginDifferenceBySupport:InvalidBounds');

    % The map output dimension must equal the ambient dimension of S.
    verifyError(testCase, @() pontryaginDifferenceBySupport( ...
        ones(2, 3), [1; 1], ones(2, 1), W), ...
        'cis2m:pontryaginDifferenceBySupport:InvalidMapRows');

    % The map input dimension must equal the ambient dimension of W.
    verifyError(testCase, @() pontryaginDifferenceBySupport( ...
        ones(2, 2), [1; 1], ones(2, 2), W), ...
        'cis2m:pontryaginDifferenceBySupport:InvalidMapColumns');
end

function testRejectsInvalidDisturbanceSets(testCase)
    G = [1; -1];
    F = [1; 1];

    % W must be one MPT Polyhedron, not another type or a set array.
    verifyError(testCase, @() pontryaginDifferenceBySupport(G, F, 1, 1), ...
        'cis2m:pontryaginDifferenceBySupport:InvalidDisturbance');
    W = Polyhedron('lb', -1, 'ub', 1);
    verifyError(testCase, @() pontryaginDifferenceBySupport(G, F, 1, [W W]), ...
        'cis2m:pontryaginDifferenceBySupport:InvalidDisturbance');

    % Support subtraction requires a nonempty compact disturbance set.
    empty_W = Polyhedron('A', [1; -1], 'b', [-1; -1]);
    verifyError(testCase, @() pontryaginDifferenceBySupport( ...
        G, F, 1, empty_W), ...
        'cis2m:pontryaginDifferenceBySupport:InvalidDisturbance');
    unbounded_W = Polyhedron('A', 1, 'b', 1);
    verifyError(testCase, @() pontryaginDifferenceBySupport( ...
        G, F, 1, unbounded_W), ...
        'cis2m:pontryaginDifferenceBySupport:InvalidDisturbance');
end

function rows = normalizedRandomRows(row_count, column_count)
    rows = randn(row_count, column_count);
    rows = rows ./ vecnorm(rows, 2, 2);
end

function [rows, bounds] = randomBoxCornerCuts(lower, upper, row_count)
    rows = normalizedRandomRows(row_count, numel(lower));

    % The support is attained at a box corner selected by each row's signs.
    % Scaling it by a factor below one guarantees that the new halfspace cuts
    % that corner, while a positive factor keeps the origin strictly feasible.
    box_support = sum(max(rows.*lower', rows.*upper'), 2);
    support_fraction = 0.55 + 0.3*rand(row_count, 1);
    bounds = support_fraction.*box_support;
end
