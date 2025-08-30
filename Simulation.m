n = 100;
x = linspace(0, 1, n)';

% Ground truth signal
u_true = sin(2 * pi * x);

% Measurement matrix and noisy observation
A = randn(50, n);
b = A * u_true + 0.01 * randn(50, 1);  % Add noise

% Reconstruct noisy signal using pseudo-inverse (just for visualization)
u_noisy = pinv(A) * b;

% Regularization setup using 1D second derivative (finite difference)
h = x(2) - x(1);
e = ones(n,1);
L = spdiags([e -2*e e], -1:1, n, n) / h^2;
L(1,:) = 0; L(end,:) = 0;  % Zero rows at boundaries

lambda = ones(n, 1) * 1e-2;        % Uniform regularization
epsilon = 1e-6;
tol = 1e-2;

% Solve the regularized problem
[u,k] = UPen(A,b,L,lambda,epsilon,tol);

% Plot the result
figure;
plot(x, u_true, 'g-', 'LineWidth', 1.5); hold on;
plot(x, u_noisy, 'r:', 'LineWidth', 1.5);
plot(x, u, 'b--', 'LineWidth', 1.5);
legend('Ground Truth', 'Noisy Signal (Pseudo-inverse)', 'Recovered Solution');
xlabel('x'); ylabel('u(x)');
title('Ground Truth, Noisy Observation, and Regularized Solution');
grid on;