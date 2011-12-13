% function [exitFlag,output] = checkDerivs(Prob, x_k, PriLev, ObjDerLev, ConsDerLev, AbsTol)
%
% TOMLAB routine for verifying user supplied routines
%
% INPUT PARAMETERS
% Prob         Problem structure created with assign routine.
% x_k          Point the check derivatives for. Default x_0 or (x_L+x_U)/2.
%              x_L and x_U have to be within 1e5.
% PriLev       Print Level, default 1. (0-1 valid).
% ObjDerLev    Depth for objective derivative check, 1 - checks gradient, 2
%              checks gradient and Hessian. Default 2 or level of derivatives 
%              supplied.
% ConsDerLev   Depth for constraint derivative check, 1 - checks Jacobian,
%              2 checks Jacobian and 2nd part of the Hessian to the
%              Lagrangian function. Default 2 or level of derivatives 
%              supplied.
% AbsTol       Absolute tolerance for errors. 
%              Default [1e-5 1e-3 1e-4 1e-3 1e-4] (g H dc d2c J).
%
% OUTPUT PARAMETERS
% exitFlag     If exitFlag ~= 0 a problem exist. See output for more
%              information. Binary indcates where problem is:
%              0 1 0 1 1. 1+2+8 = 13. Problems everywhere but 'dc', 'J'.
%              1 1 1 1 1 = 'J' 'd2c' 'dc' 'H' 'g'.
% output       Structure containing analysis information.
%   g,H,dc,d2c,J   Structure with results.
%     minErr     The smallest error.
%     avgErr     The average error.
%     maxErr     The largest error.
%     idx        Index for elements with errors.
%     exitFlag   1 if problem with the function.

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2006 by Tomlab Optimization Inc., $Release: 5.7.0$
% Written Aug 10, 2004.   Last modified Dec 26, 2006.

function [exitFlag,output] = checkDerivs(Prob, x_k, PriLev, ObjDerLev, ConsDerLev, AbsTol)

if nargin < 1
    error('checkDerivs requires at least one parameter, Prob'); 
end

if nargin < 6
    AbsTol=[];
    if nargin < 5
        ConsDerLev=2;
        if nargin < 4
            ObjDerLev=2;
            if nargin < 3
                PriLev = 1;
                if nargin < 2
                    x_k = [];
                end
            end
        end
    end
end

if isempty(PriLev)
    PriLev = 1;
end

if isempty(ObjDerLev)
    ObjDerLev = 2;
end

if isempty(ObjDerLev)
    ConsDerLev = 2;
end

if isempty(AbsTol)
    AbsTol=[1e-5; 1e-3; 1e-4; 1e-3; 1e-4];
end

if length(AbsTol) ~= 5
    error('AbsTol must be of length 5');
end

% Start with check on Patterns, remove LargeScale if set.

if isfield(Prob, 'LargeScale')
    Prob.LargeScale = 0; % Makes sure iniSolve doesn't estimate patterns
end

% Check the Prob structure

Prob = ProbCheck(Prob, 'conSolve'); % Default solver.
Prob = iniSolve(Prob,12,ObjDerLev,ConsDerLev);

% Check the patterns before anything else, error is returned if incorrect

if ~isempty(Prob.ConsPattern)
    ConsPattern = estConsPattern(Prob);
    if ~isempty(find(sparse(ConsPattern - Prob.ConsPattern) ~= 0))
        error('ConsPattern is incorrect, cannot proceed');
    end
end

if ~isempty(Prob.JacPattern)
    JacPattern = estJacPattern(Prob);
    if ~isempty(find(sparse(JacPattern - Prob.JacPattern) ~= 0))
        error('JacPattern is incorrect, cannot proceed');
    end
end

if ~isempty(Prob.HessPattern)
    HessPattern = estHessPattern(Prob);
    if ~isempty(find(sparse(HessPattern - Prob.HessPattern) ~= 0))
        error('HessPattern is incorrect, cannot proceed');
    end
end

% Generate output
oF = struct('exitFlag', 0, 'minErr', [], 'avgErr', [], 'maxErr', [], 'idx', []);
output = struct('g', oF, 'H', oF, 'dc', oF, 'd2c', oF, 'J', oF);
clear oF;
exitFlag = 0;      % Default OK.

if isfield(Prob, 'PartSep')
    if isfield(Prob.PartSep, 'pSepFunc')
        if Prob.PartSep.pSepFunc == 2
            fprintf('checkDerivs cannot be used for partially separable functions. \n');
            return;
        end    
    end
end

% Setup ObjDerLev and ConsDerLev, check if least squares problem.

FUNCS = Prob.FUNCS;

objSwitch  = 0; % 0 standard, 1 sim problem, 2 least squares.

if ~isempty(FUNCS.fc)
    objSwitch  = 1; % A sim problem.
elseif ~isempty(FUNCS.r)
    objSwitch = 2; % A least squares problem.
end

% Switching ObjDerLev and ConsDerLev depending on what is given in
% Prob.FUNCS

% Standard case

if objSwitch == 0
    if isempty(FUNCS.g) & isempty(FUNCS.H)
        ObjDerLev = 0;
    elseif isempty(FUNCS.H)
        ObjDerLev = 1;
    end
    if isempty(FUNCS.dc) & isempty(FUNCS.d2c)
        ConsDerLev = 0;
    elseif isempty(FUNCS.d2c)
        ConsDerLev = 1;
    end
end

% sim problems only 1 level

if objSwitch == 1
    if isempty(FUNCS.gdc)
        ObjDerLev  = 0;
        ConsDerLev = 0;        
    else
        ObjDerLev  = 1;
        ConsDerLev = 1;        
    end
end

if objSwitch == 2
    if isempty(FUNCS.J) & isempty(FUNCS.d2r)
        ObjDerLev = 0;
    elseif isempty(FUNCS.d2r)
        ObjDerLev = 1;
    end
    if isempty(FUNCS.dc) & isempty(FUNCS.d2c)
        ConsDerLev = 0;
    elseif isempty(FUNCS.d2c)
        ConsDerLev = 1;
    end
end


if ObjDerLev == 0 & ConsDerLev == 0 & PriLev > 0
    fprintf('No derivatives supplied, no checks done. \n');
    return;
end

% Check point for checking derivatives.

if isempty(Prob.x_L)
    Prob.x_L = -inf*ones(Prob.N,1);
end

if isempty(Prob.x_U)
    Prob.x_U = inf*ones(Prob.N,1);
end

if isempty(x_k)
    if ~isempty(Prob.x_0)
        x_k = Prob.x_0(:);
    else
        x_L = max(Prob.x_L(:), -1e5);
        x_U = min(Prob.x_U(:),  1e5);
        x_k = (x_L + x_U) / 2;
    end
end

if objSwitch == 0
    if ObjDerLev >= 1
        % Check gradient first, if OK use for Hessian check.
        Prob.NumDiff = [];
        fAnal = nlp_f(x_k, Prob);
        gAnal = nlp_g(x_k, Prob);
        Prob.NumDiff = 1;
        global NLP_x NLP_f NLP_g NLP_xg
        NLP_x=[]; NLP_f=[]; NLP_g=[]; NLP_xg=[];
        gNum  = nlp_g(x_k, Prob);
        if ~all(size(gAnal) == size(gNum))
            fprintf('Analytical gradient not same size as numerical, cannot proceed. \n');
            return
        end        
        gErr = abs(gAnal - gNum);
        output.g.minErr = min(gErr);
        output.g.maxErr = max(gErr);
        if ~isempty(gErr)
            output.g.avgErr = mean(gErr);
        end
        output.g.idx    = find(gErr>AbsTol(1));
        if output.g.maxErr > AbsTol(1)
            exitFlag = exitFlag + 1;
            output.g.exitFlag = 1;
        end
        if PriLev > 0
            fprintf('Largest error in gradient is %0.6g.',full(output.g.maxErr));
            fprintf('\n');            
        end
    end
    if ObjDerLev == 2 % Check gradient and Hessian.
        if output.g.exitFlag == 0 % Use gradient for Hessian
            NumDiff = -1;
        else
            NumDiff =  1;
        end
        Prob.NumDiff = [];
        Prob.HessPattern = triu(Prob.HessPattern);
        HAnal = nlp_H(x_k, Prob);
        Prob.NumDiff = NumDiff;            
        global NLP_g NLP_xg NLP_pSepIndex NLP_x NLP_f
        NLP_g=[]; NLP_xg=[]; NLP_pSepIndex =[]; NLP_x=[]; NLP_f=[];
        HNum  = nlp_H(x_k, Prob);
        if ~all(size(HAnal) == size(HNum))
            fprintf('Analytical Hessian not same size as numerical, cannot proceed. \n');
            return
        end
        HErr = abs(HAnal - HNum);
        output.H.minErr = min(min(HErr));
        output.H.maxErr = max(max(HErr));
        if ~isempty(HErr)
            output.H.avgErr = mean(mean(HErr));
        end
        output.H.idx = find(HErr>AbsTol(2));        
        if output.H.maxErr > AbsTol(2)
            exitFlag = exitFlag + 2;
            output.H.exitFlag = 1;
        end        
        if PriLev > 0
            fprintf('Largest error in Hessian is %0.6g.',full(output.H.maxErr));
            fprintf('\n');
        end
    end
end

if any([0,2] == objSwitch)
    if ConsDerLev >= 1
        % Check Jacobian first, if OK use for d2c check.
        Prob.ConsDiff = [];
        dcAnal = nlp_dc(x_k, Prob);
        global NLP_xc NLP_xdc NLP_c NLP_dc
        NLP_xc=[]; NLP_xdc=[]; NLP_c=[]; NLP_dc=[];
        Prob.ConsDiff = 1;
        dcNum  = nlp_dc(x_k, Prob);
        if ~all(size(dcAnal) == size(dcNum))
            fprintf('Analytical constraint Jacobian not same size as numerical, cannot proceed. \n');
            return            
        end 
        dcErr = abs(dcAnal - dcNum);
        output.dc.minErr = min(min(dcErr));
        output.dc.maxErr = max(max(dcErr));
        if ~isempty(dcErr)
            output.dc.avgErr = mean(mean(dcErr));
        end
        output.dc.idx    = find(dcErr>AbsTol(3));
        if output.dc.maxErr > AbsTol(3)
            exitFlag = exitFlag + 4;
            output.dc.exitFlag = 1;
        end
        if PriLev > 0
            fprintf('Largest error in constraint Jacobian is %0.6g.',full(output.dc.maxErr));
            fprintf('\n');
        end        
    end
    if ConsDerLev == 2 % Check d2c Prob.mNonLin used
        if output.dc.exitFlag == 0 % Use Jacobian for d2c check
            ConsDiff = -1;
        else
            ConsDiff =  1;
        end
        Prob.ConsDiff = [];
        lam = rand(Prob.mNonLin,1);
        d2cAnal = nlp_d2c(x_k, lam, Prob);
        Prob.ConsDiff = ConsDiff;
        global NLP_xdc
        NLP_xdc=[];
        global NLP_xc NLP_xdc NLP_c NLP_dc
        NLP_xc=[]; NLP_xdc=[]; NLP_c=[]; NLP_dc=[];
        d2cNum  = nlp_d2c(x_k, lam, Prob);
        if ~all(size(d2cAnal) == size(d2cNum))
            fprintf('Analytical d2c not same size as numerical, cannot proceed. \n');
            return            
        end
        d2cErr = abs(d2cAnal - d2cNum);
        output.d2c.minErr = min(min(d2cErr));
        output.d2c.maxErr = max(max(d2cErr));
        if ~isempty(d2cErr)
            output.d2c.avgErr = mean(mean(d2cErr));
        end
        output.d2c.idx    = find(d2cErr>AbsTol(4));        
        if output.d2c.maxErr > AbsTol(4)
            exitFlag = exitFlag + 8;
            output.d2c.exitFlag = 1;
        end
        if PriLev > 0
            fprintf('Largest error in d2c is %0.6g.',full(output.d2c.maxErr));
            fprintf('\n');
        end               
    end
    
    if objSwitch == 2 % Check J
        if ObjDerLev >= 1
            % Check Jacobian first, if OK use for d2r check.
            Prob.NumDiff = [];
            JAnal = nlp_J(x_k, Prob);
            global LS_x LS_r LS_xJ LS_J
            LS_x=[]; LS_r=[]; LS_xJ=[]; LS_J=[];
            Prob.NumDiff = 1;
            JNum  = nlp_J(x_k, Prob);
            if ~all(size(JAnal) == size(JNum))
                fprintf('Analytical residual Jacobian not same size as numerical, cannot proceed. \n');
                return
            end
            JErr = abs(JAnal - JNum);
            output.J.minErr = min(min(JErr));
            output.J.maxErr = max(max(JErr));
            if ~isempty(JErr)
                output.J.avgErr = mean(mean(JErr));
            end
            output.J.idx    = find(JErr>AbsTol(5));
            if output.J.maxErr > AbsTol(5)
                exitFlag = exitFlag + 16;
                output.J.exitFlag = 1;
            end
            if PriLev > 0
                fprintf('Largest error in residual Jacobian is %0.6g.',full(output.J.maxErr));
                fprintf('\n');
            end                        
        end
    end
end

% MODIFICATION LOG
%
% 040811 med   Created.
% 040811 med   Fixed defaults.
% 040811 med   Default for AbsTol corrected.
% 040811 med   Added printing of results.
% 040811 med   Error if not same sized on analytical and numerical results.
% 040811 med   fprintf not handling sparse inputs.
% 040811 med   Error in residual Jacobian check.
% 040811 med   Partially separate functions removed.
% 040922 med   Part separable function isfield check added.
% 041129 med   Test for Patterns added.
% 050902 med   d2c check updated (rand for lam's now)
% 060814 med   FUNCS used for callbacks instead
% 061226 med   HessPattern corrected