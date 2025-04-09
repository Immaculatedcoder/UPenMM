close all
clear
clc
%%%
n = 100;
kappa = 1;
[A,b,u] = heat(n,kappa);

%%% noisy data %%%
eta = rand(n,1);
eta = eta/norm(eta);

delta = 0.1; % Noise level

noise = delta * eta * norm(A*b);
b_noise = b + noise;

figure(1);
plot(b)
hold on
plot(b_noise)
legend('Original data', 'Noisy data')
hold off




% What we have so far
% A = Matrix operator
% b_noise = given sigal
% u = Signal to be recovered.

L = sdo(n);
epsilon = 1e-6;
tol = 1e-5;

lambda = zeros(n,1);

[u_p,k] = UPen(A,b,L,lambda,epsilon,tol);
disp(k)
figure(2)
plot(u,'r-'); % What we wish to recover.
hold on
plot(u_p(1:end-1),'k--',LineWidth=2)
legend('exact','UPen')
hold off






