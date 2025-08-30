% %%%%%%%%%%%%%%%%%%%%%%%
n = 100;
x = linspace(0, 1, n+1); x(end) = [];
y = sin(2*pi*x);                          % signal
y_dd_exact = -4*pi^2*sin(2*pi*x);         % true second derivative

D2 = sdo(n);
y_dd = (D2 * y')';                        % numerical second derivative

plot(x, y_dd, 'r--', x, y_dd_exact, 'b-');
legend('Numerical d^2y/dx^2', 'Exact d^2y/dx^2');
title('Second Derivative via Finite Difference');
xlabel('x'); ylabel('d^2y/dx^2');
grid on

% Error
max_error = max(abs(y_dd - y_dd_exact));
disp(['Max error: ', num2str(max_error)])




% 
% function [u, obj_val] = solve_regularized_least_squares(A, b, L, lambda, epsilon)
% % Solves: min_u (1/2)||Au - b||^2 + sum_i lambda_i * [(Lu)_i^2 + epsilon]
% %
% % Inputs:
% %   A      - m x n matrix
% %   b      - m x 1 vector
% %   L      - N x n matrix (e.g., finite difference second derivative)
% %   lambda - N x 1 vector of regularization weights
% %   epsilon- small constant (scalar, e.g. 1e-6)
% %
% % Outputs:
% %   u        - solution vector (n x 1)
% %   obj_val  - value of the objective function at solution u
% 
%     % Compute the regularization matrix
%     Lambda = diag(lambda);                % N x N
%     Reg = L' * Lambda * L;                % n x n (symmetric PSD)
% 
%     % Construct the normal equation matrix
%     H = A' * A + 2 * Reg;                 % Hessian-like term
%     rhs = A' * b;                         % Right-hand side
% 
%     % Solve the linear system
%     u = H \ rhs;                          % Could also use pcg() for large problems
% 
%     % Compute the objective function value including epsilon
%     Lu = L * u;
%     obj_val = 0.5 * norm(A * u - b)^2 + sum(lambda .* (Lu.^2 + epsilon));
% end
% 
% 
% n = 100;
% x = linspace(0, 1, n)';
% 
% % Ground truth signal
% u_true = sin(2 * pi * x);
% 
% % Measurement matrix and noisy observation
% A = randn(50, n);
% b = A * u_true + 0.01 * randn(50, 1);  % Add noise
% 
% % Reconstruct noisy signal using pseudo-inverse (just for visualization)
% u_noisy = pinv(A) * b;
% 
% % Regularization setup using 1D second derivative (finite difference)
% h = x(2) - x(1);
% e = ones(n,1);
% L = spdiags([e -2*e e], -1:1, n, n) / h^2;
% L(1,:) = 0; L(end,:) = 0;  % Zero rows at boundaries
% 
% lambda = ones(n, 1) * 1e-2;        % Uniform regularization
% epsilon = 1e-6;
% 
% % Solve the regularized problem
% [u, obj_val] = solve_regularized_least_squares(A, b, L, lambda, epsilon);
% 
% % Plot the result
% figure;
% plot(x, u_true, 'g-', 'LineWidth', 1.5); hold on;
% plot(x, u_noisy, 'r:', 'LineWidth', 1.5);
% plot(x, u, 'b--', 'LineWidth', 1.5);
% legend('Ground Truth', 'Noisy Signal (Pseudo-inverse)', 'Recovered Solution');
% xlabel('x'); ylabel('u(x)');
% title('Ground Truth, Noisy Observation, and Regularized Solution');
% grid on;
% 
