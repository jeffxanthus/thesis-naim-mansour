% function r = ls1_r(x, Prob)
%
% Compute residuals to nonlinear least squares problem Gisela

function r = ls1_r(x, Prob)

% This is the standard TOMLAB global parameter to be used in the
% communication between the residual and the Jacobian routine

global US_A US_B

% The extra weight parameter K is sent as part of the structure
K  = Prob.userParam.K;

t  = Prob.LS.t(:);   % Pick up the time points

% Exponential computations takes time, and may be done once, and
% reused when computing the Jacobian
US_A = exp(-x(1)*t);
US_B = exp(-x(2)*t);

% The loop version is used in ls_prob.m, problem 3
% It is more efficient to compute without looping:

r = K*x(1)*(US_B - US_A) / (x(3)*(x(1)-x(2))) - Prob.LS.y;

% This is the code to compute the residual written as a loop
% m  = size(Prob.LS.y,1);
% for i=1:m
%     r(i)=K*x(1)*(US_B(i)-US_A(i))/(x(3)*(x(1)-x(2)));
% end
% r = r(:) - Prob.LS.y;