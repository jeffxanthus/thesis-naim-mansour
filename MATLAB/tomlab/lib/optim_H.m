% optim_H.m
%
% function H = optim_H(x, Prob, varargin)
%
% optim_H is used to implement the OPT TB 2.0 interface
%
% The Hessian matrix H is returned, if available in the global variable NLP_H
% with the corresponding x value in NLP_xH
%
% optim_H is called from the TOMLAB gateway function nlp_H.
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written July 29, 1999.     Last modified Jan 15, 2003.

function H = optim_H(x, Prob)

global NLP_xg NLP_g NLP_xH NLP_H

% Return the Hessian, if computed.

if ~isempty(NLP_xH)
    if all(x==NLP_xH)
        H=NLP_H;
    else
        f=optim_fgH(x, Prob);
        H=NLP_H;
    end
else
    H=[];
end

% MODIFICATION LOG:
%
% 030115 hkh Clean routine after major revision of fmincon
