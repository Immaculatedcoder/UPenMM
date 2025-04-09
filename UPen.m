function [u,k] = UPen(A,b,lambda,tol)
%-------------------------------------------------------------------------%
% Uniform Penalty Based on Majorization-Minimization Method
% 
% INPUTS:
%   A   - Matrix
%   b   - Noisy data
%   lambda   - legularization vector of dimension p 
%
% OUTPUT:
%   u      - Solution
%   k    - Number of iterations
% USAGE:
%   [coeffs, f_haar_c, f_haar_x] 
%-------------------------------------------------------------------------%
% Author: Emmanuel Adebayo
% Email: adebayo@udel.edu
% Date: 08-April-2025
%-------------------------------------------------------------------------%


k=0;
p = length(lambda);
while k < 5
    pre_lambda = lambda;
    u = preglsq(A,b,lambda);
    
    for i = 1:p
        lambda(i) = (norm(A*u - b)^2)/(p*psi(i));
    end

    % Checking for Convergence
    if norm(lambda - pre_lambda) <= norm(pre_lambda)*tol
        break;
    end
    
    k = k + 1;
end
end