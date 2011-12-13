% optim_g.m
%
% function g = optim_g(x, Prob, varargin)
%
% optim_g is used to implement the OPT TB 2.0 interface
%
% The gradient g is returned, if available in the global variable NLP_g
% with the corresponding x value in NLP_xg
%
% optim_g is called from the TOMLAB gateway function nlp_g.

% Kenneth Holmstrom, Tomlab Optimization AB, E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2005 by Tomlab Optimization AB, Sweden. $Release: 4.7.0$
% Written July 29, 1999.     Last modified Jan 15, 2003.

% function g = optim_g(x, Prob, varargin)

function g = optim_g(x, Prob)

global NLP_xg NLP_g NLP_xH NLP_H

% Return the gradient, if computed.

if ~isempty(NLP_xg)
    if all(x==NLP_xg)
        g=NLP_g;
    else
        f=optim_fgH(x, Prob);
        g=NLP_g;
    end
else
    g=[];
end

% MODIFICATION LOG:
%
% 030115 hkh Clean routine after major revision of fmincon
