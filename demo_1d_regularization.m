function demo_1d_regularization()
    % Parameters
    N = 100;                % Number of grid points
    x = linspace(0, 1, N)'; % Domain
    h = x(2) - x(1);        % Grid spacing

    % Ground truth signal
    u_true = sin(2*pi*x) + 0.1 * randn(N,1);  % Noisy sinusoidal signal

    % Measurement matrix A and observation b
    A = eye(N);             % Identity = observe all points
    b = u_true + 0.05*randn(N, 1);  % Add noise

    % Regularization matrix (2nd derivative operator)
    L = build_L_1d(N, h);   % Finite difference 2nd derivative
    lambda = 1e-2 * ones(N,1);  % Regularization weights
    epsilon = 1e-6;             % Stability constant (can be ignored in quad form)

    % Solve the regularized least squares problem
    u_sol = solve_regularized_least_squares(A, b, L, lambda, epsilon);

    % Plot
    figure;
    plot(x, u_true, 'g-', 'LineWidth', 1.5); hold on;
    plot(x, b, 'rx', 'DisplayName', 'Noisy data');
    plot(x, u_sol, 'b-', 'LineWidth', 2, 'DisplayName', 'Smoothed solution');
    legend('True signal', 'Noisy data', 'Recovered u');
    xlabel('x'); ylabel('u(x)');
    title('1D Regularized Least Squares with 2nd Derivative Penalty');
end

function L = build_L_1d(N, h)
    e = ones(N,1);
    L = spdiags([e -2*e e], -1:1, N, N);
    L(1,:) = 0; L(end,:) = 0;  % Handle boundaries (zero for now)
    L = L / h^2;
end

function u = solve_regularized_least_squares(A, b, L, lambda, epsilon)
    Lambda = diag(lambda);         % Diagonal weighting
    Reg = L' * Lambda * L;         % Regularization matrix
    H = A' * A + 2 * Reg;          % Composite system matrix
    rhs = A' * b;                  % Right-hand side
    u = H \ rhs;                   % Solve system
end
