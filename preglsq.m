function u = preglsq(A, b, L, lambda)
%PREGLSQ Summary of this function goes here
%   Detailed explanation goes here
    Lambda = diag(lambda);                % N x N
    Reg = L' * Lambda * L;                % n x n (symmetric PSD)

    % Construct the normal equation matrix
    H = A' * A + 2 * Reg;                 % Hessian-like term
    rhs = A' * b;                         % Right-hand side

    % Solve the linear system
    u = H \ rhs;
end

