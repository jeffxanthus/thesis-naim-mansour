% TOMLAB - Financial Optimization, Mean Variance. Example 1

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004 by Tomlab Optimization Inc., Sweden. $Release: 4.4.0$
% Written Oct 16, 2004. Last modified Oct 16, 2004.
%

% We have 50 assets.
% We have the returns for 60 periods.

% The first 25 assets have been going up with 0.001-0.025 % per period
% plus some random number.

clear all

R = zeros(60,50);

rand('state',0);

for i=1:25
    R(1:60,i) = 0.001*i*ones(60,1)+rand(60,1)/1000;
end

% The last 25 assets have been going down with 0.001-0.025 % per period

for i=26:50
    R(1:60,i) = -0.001*i*ones(60,1)+rand(60,1)/1000;
end

r_min = 0.015; % MINIMUM RETURN

w_min = zeros(50,1);
w_max = 0.5*ones(50,1);

solver = 'snopt';

% TESTING MODEL 1.

w1 = portoptMeanVariance(R, w_min, w_max, r_min, solver, [], 1, 'MYOPT1', 1);

% TESTING MODEL 2.

w2 = portoptMeanVariance(R, w_min, w_max, r_min, solver, [], 1, 'MYOPT2', 2);

% Moving down minimum return to 0.005

r_min = 0.005;
w3 = portoptMeanVariance(R, w_min, w_max, r_min, solver, [], 1, 'MYOPT3', 2);

% Moving up minimum return to 0.022

r_min = 0.022;
w4 = portoptMeanVariance(R, w_min, w_max, r_min, solver, [], 1, 'MYOPT4', 2);

% We must invest at least 1.5 % in each asset positive asset, r_min = 0.003

r_min = 0.010;
w_min(1:25,1) = 0.015*ones(25,1);
w5 = portoptMeanVariance(R, w_min, w_max, r_min, solver, [], 1, 'MYOPT5', 2);

% We cannot invest anything in first 5 assets.

r_min = 0.010;
w_min(1:25,1) = 0.015*ones(25,1);
w_min(1:5,1) = 0*ones(5,1);
w_max(1:5,1) = w_min(1:5,1);
w6 = portoptMeanVariance(R, w_min, w_max, r_min, solver, [], 1, 'MYOPT6', 2);

format short
[w1 w2 w3 w4 w5 w6 [1:50]']

% MODIFICATION LOG
%
% 041016  med  Created