function tests = test_transform_to_brunovsky_normal_form
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    matlab_dir = fileparts(fileparts(mfilename('fullpath')));
    testCase.addTeardown(@path, path);
    addpath(matlab_dir, fullfile(matlab_dir, 'support_functions'));

    % The transformation uses ctrb from Control System Toolbox.
    assumeTrue(testCase, exist('ctrb', 'file') == 2, ...
        'Requires Control System Toolbox.');
end

function testFullyActuatedSystemPreservesDynamicsAndConstraints(testCase)
    % A nonsingular square B gives two controllability chains of length one.
    A = [2 1; 0 3];
    B = [2 1; 1 1];
    E = [1; -2];
    Gxu = [eye(4); -eye(4); 1 -2 3 -4];
    Fxu = (1:9)';
    [Ac, Bc, Ec, Gc, Fc, T, nu, Am, Bm] = ...
        transformToBrunovskyNormalForm(A, B, E, Gxu, Fxu);

    % Fully actuated Brunovsky dynamics have Ac = 0 and nilpotency index one.
    verifySize(testCase, Ac, [2 2]);
    verifySize(testCase, Bc, [2 2]);
    verifySize(testCase, Gc, [9 4]);
    verifyEqual(testCase, nu, 1);
    verifyEqual(testCase, Ac, zeros(2), 'AbsTol', 1e-12);
    verifyEqual(testCase, Fc, Fxu);

    % Verify all transformation identities over three state-input-disturbance
    % samples instead of checking only the returned matrices independently.
    x = [1 -2 3; 4 2 -1];
    u = [-3 1 2; 2 -1 4];
    w = [0.5 -1 2];
    verifyTransformationIdentities(testCase, A, B, E, Gxu, ...
        Ac, Bc, Ec, Gc, T, Am, Bm, x, u, w, 1e-11);
end

function testUnequalChainLengthsRetainOriginalNilpotency(testCase)
    % This system has controllability indices two and one.
    A = [0 1 0; 1 2 0; 0 0 3];
    B = [0 0; 1 0; 0 1];
    Gxu = [eye(5); -eye(5)];
    [Ac, Bc, Ec, Gc, ~, ~, nu] = ...
        transformToBrunovskyNormalForm(A, B, [], Gxu, ones(10, 1));

    verifySize(testCase, Ac, [3 3]);
    verifySize(testCase, Bc, [3 2]);
    verifySize(testCase, Gc, [10 5]);

    % An absent disturbance map remains absent after the coordinate change.
    verifyEmpty(testCase, Ec);

    % Ac^nu = 0 and Ac^(nu-1) ~= 0 establish that nu is the exact index.
    verifyEqual(testCase, nu, 2);
    verifyEqual(testCase, Ac^nu, zeros(3), 'AbsTol', 1e-12);
    verifyGreaterThan(testCase, norm(Ac^(nu - 1), 'fro'), 0);
end

function testSparseInputsAndEmptyConstraints(testCase)
    % computeRCIS calls this function with sparse data. A 0-by-(n+m) Gxu and
    % 0-by-1 Fxu are the canonical representation of no joint constraints.
    A = sparse([0 1; -2 3]);
    B = sparse([0; 1]);
    Gxu = sparse([], [], [], 0, 3);
    Fxu = sparse([], [], [], 0, 1);

    [Ac, Bc, Ec, Gc, Fc, T, nu, Am, Bm] = ...
        transformToBrunovskyNormalForm(A, B, [], Gxu, Fxu);

    verifyEqual(testCase, nu, 2);
    verifyTrue(testCase, issparse(Ac));
    verifyTrue(testCase, issparse(Bc));
    verifyEqual(testCase, full(Ac), [0 1; 0 0], 'AbsTol', 1e-12);
    verifyEqual(testCase, full(Bc), [0; 1], 'AbsTol', 1e-12);
    verifyEmpty(testCase, Ec);
    verifySize(testCase, Gc, [0 3]);
    verifySize(testCase, Fc, [0 1]);

    % With no disturbance, the nominal dynamics identity still must hold.
    x = [1 -2 3; 4 2 -1];
    u = [-3 1 2];
    z = T*x;
    r = Am*z + Bm*u;
    verifyEqual(testCase, rank(full(T)), size(T, 1));
    verifyEqual(testCase, rank(full(Bm)), size(Bm, 1));
    verifyEqual(testCase, T \ z, x, 'AbsTol', 1e-12);
    verifyEqual(testCase, Bm \ (r - Am*z), u, 'AbsTol', 1e-12);
    verifyEqual(testCase, Ac*z + Bc*r, T*(A*x + B*u), 'AbsTol', 1e-12);
end

function testWarnsForIllConditionedControllabilityBasis(testCase)
    % The system remains controllable, but differently scaled input channels
    % make the controllability basis exceed the diagnostic threshold.
    A = zeros(2);
    B = diag([1 1e-15]);
    Gxu = [eye(4); -eye(4)];
    Fxu = ones(8, 1);

    lastwarn('');
    [Ac, Bc, ~, ~, ~, T, nu, ~, Bm] = ...
        transformToBrunovskyNormalForm(A, B, [], Gxu, Fxu);
    [warning_message, warning_id] = lastwarn;

    % The warning currently has no identifier, so verify its exact message.
    verifyEmpty(testCase, warning_id);
    verifyEqual(testCase, warning_message, 'Condition number > 1e14.');

    % Ill-conditioning must not alter the canonical result for this system.
    verifyEqual(testCase, nu, 1);
    verifyEqual(testCase, Ac, zeros(2), 'AbsTol', 1e-12);
    verifyEqual(testCase, Bc, eye(2), 'AbsTol', 1e-12);
    verifyEqual(testCase, rank(full(T)), 2);
    verifyEqual(testCase, rank(full(Bm)), 2);
end

function testRejectsUncontrollableSystem(testCase)
    % The third state is unreachable, so the controllability matrix has rank
    % two while the state dimension is three.
    A = diag([0 1 2]);
    B = [1; 1; 0];
    Gxu = [eye(4); -eye(4)];
    Fxu = ones(8, 1);

    verifyError(testCase, ...
        @() transformToBrunovskyNormalForm(A, B, [], Gxu, Fxu), ...
        'cis2m:transformToBrunovskyNormalForm:UncontrollableSystem');
end

function testSevenStateThreeInputSystem(testCase)
    % Build an original-coordinate system from three Brunovsky chains with
    % controllability indices 4, 2, and 1. The nontrivial similarity transform
    % and feedback make every returned identity meaningful.
    chain_lengths = [4 2 1];
    [Ac_expected, Bc_expected] = brunovskyChains(chain_lengths);
    n = size(Ac_expected, 1);
    m = size(Bc_expected, 2);

    T_seed = tril(ones(n)) + diag(1:n);
    Bm_seed = diag([2 -1 3]);
    Am_seed = [1 -2 0 1 0 1 -1;
               0  1 2 0 1 0  1;
              -1  0 1 2 0 1  3];
    Ec_seed = [1 0; 0 -1; 2 1; -1 2; 1 1; 0 2; -2 1];
    A = T_seed \ ((Ac_expected + Bc_expected*Am_seed)*T_seed);
    B = T_seed \ (Bc_expected*Bm_seed);
    E = T_seed \ Ec_seed;

    Gxu = [eye(n + m); -eye(n + m);
        1 -2 0 3 0 -1 2 4 -3 1;
        0 1 -1 0 2 3 -2 -1 2 4;
        2 0 1 -2 1 0 3 2 1 -3];
    Fxu = (10:(9 + size(Gxu, 1)))';
    [Ac, Bc, Ec, Gc, Fc, T, nu, Am, Bm] = ...
        transformToBrunovskyNormalForm(A, B, E, Gxu, Fxu);

    % Recover the prescribed canonical structure and full-rank transforms.
    verifyEqual(testCase, nu, 4);
    verifySize(testCase, Ac, [n n]);
    verifySize(testCase, Bc, [n m]);
    verifySize(testCase, Gc, size(Gxu));
    verifyEqual(testCase, rank(full(T)), n);
    verifyEqual(testCase, rank(full(Bm)), m);
    verifyEqual(testCase, Ac, Ac_expected, 'AbsTol', 5e-10);
    verifyEqual(testCase, Bc, Bc_expected, 'AbsTol', 5e-10);
    verifyEqual(testCase, Ac^nu, zeros(n), 'AbsTol', 5e-10);
    verifyGreaterThan(testCase, norm(Ac^(nu - 1), 'fro'), 0.5);
    verifyEqual(testCase, Fc, Fxu);

    % Check the identities over four nontrivial samples in seven dimensions.
    x = reshape(sin(1:(4*n)), n, 4);
    u = reshape(cos(1:(4*m)), m, 4);
    w = [0.2 -0.5 1 2; -1 0.25 0.5 -0.75];
    verifyTransformationIdentities(testCase, A, B, E, Gxu, ...
        Ac, Bc, Ec, Gc, T, Am, Bm, x, u, w, 2e-9);
end

function verifyTransformationIdentities(testCase, A, B, E, Gxu, ...
        Ac, Bc, Ec, Gc, T, Am, Bm, x, u, w, tolerance)
    % Apply z = T*x and r = Am*z + Bm*u.
    z = T*x;
    r = Am*z + Bm*u;

    % T and Bm must be invertible for (x, u) <-> (z, r) to be bijective.
    verifyEqual(testCase, rank(full(T)), size(T, 1));
    verifyEqual(testCase, rank(full(Bm)), size(Bm, 1));
    verifyEqual(testCase, T \ z, x, 'AbsTol', tolerance);

    % The feedback map must recover the original physical inputs.
    verifyEqual(testCase, Bm \ (r - Am*z), u, 'AbsTol', tolerance);

    % Joint constraints must evaluate identically in both coordinates.
    verifyEqual(testCase, Gc*[z; r], Gxu*[x; u], 'AbsTol', tolerance);

    % Transformed and original dynamics must describe the same next state.
    verifyEqual(testCase, Ac*z + Bc*r + Ec*w, ...
        T*(A*x + B*u + E*w), 'AbsTol', tolerance);
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
