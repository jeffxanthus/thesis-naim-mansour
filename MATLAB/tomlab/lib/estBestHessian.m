% function [g_k, H_k] = estBestHessian(Result)
%
% estBestHessian estimates the best Hessian. Result.x_k(:,1) will be used
% for the estimation. The best step-size is estimated by TOMLAB.
%
% If the gradient is given it will be used. The analytical hessian is
% returned if given.
%
% INPUT:
% Result      The Tomlab result structure, normally returned by tomRun
%
% OUTPUT:
% g_k         The gradient at Result.x_k(:,1)
% H_k         Hessian of the objective at Result.x_k(:,1)

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written Jan 31, 2005.  Last modified Jun 9, 2007.

function [g_k,H_k] = estBestHessian(Result)

if nargin < 1
    error('estBestHessian needs the TOMLAB result as input');
end

% Reset global parameters
global NLP_g NLP_xg NLP_H NLP_xH
NLP_g = [];
NLP_xg = [];
NLP_H = [];
NLP_xH = [];
global gTol HTol
% Set empty, forcing estimate of best step size
gTol = [];
HTol = [];
% Set negative, otherwise Prob.optParam.DiffInt is used instead of the
% algorithm estimating best possible step size
Result.Prob.optParam.DiffInt = -1;
Result.Prob.GradTolg = [];
Result.Prob.GradTolH = [];
Result.Prob.NumDiff = 1;
Result.Prob.NumDiff = 1;
if Result.Prob.ADObj ~= 0;
    nlp_f(Result.x_k(:,1),Result.Prob);
end
g_k = nlp_g(Result.x_k(:,1),Result.Prob);
H_k = nlp_H(Result.x_k(:,1),Result.Prob);   % The Hessian of the objective function

% MODIFICATION LOG
%
% 050131  hkh  Algorithm formulated and written
% 070609  med  Added call to nlp_f in case of MAD