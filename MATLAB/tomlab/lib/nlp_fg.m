% nlp_fg.m
%
% function [Mode, f, g]=nlp_fg(x, Prob, Mode, nState)
%
% nlp_fg calls the TOMLAB gateway function nlp_f,
% which computes the function value f(x) 
%
% nlp_fg also calls the TOMLAB gateway function nlp_g,
% which computes the gradient g at x 
%
% Mode and nState is sent as Prob.Mode, Prob.nState to nlp_f and nlp_g.
%
% Mode = 0 Assign function values
% Mode = 1 Assign known derivatives, unknown set as -11111 (=missing value)
% Mode = 2 Assign function and known derivatives, unknown derivatives -11111
%
% nState = 1         First call
% nState = 0         Other calls
% nState = 2+Inform  Last call, see Inform parameter for the solver
%
% nlp_fg is used when calling MINOS and SNOPT.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1997-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Mar 5, 1997.    Last modified July 24, 2011.

function [Mode, f, g]=nlp_fg(x, Prob, Mode, nState)

if ~isempty(Prob.FUNCS.fg)
    if Mode==0
        [Mode,f]   = feval(Prob.FUNCS.fg,x,Prob,Mode,nState);
    else
        [Mode,f,g] = feval(Prob.FUNCS.fg,x,Prob,Mode,nState);
    end

else
    % Mode and nState set into the Prob structure
    Prob.Mode   = Mode;
    Prob.nState = nState;
    f           = nlp_f(x, Prob );
    if Mode > 0
        if Prob.NumDiff ~= 6
            g   = nlp_g(x, Prob );
        else
            g   = [];
        end
    end
end

% MODIFICATION LOG:
%
% 990626  hkh  Avoid feval
% 000820  hkh  Revised for v3.0, changing arguments
% 000828  hkh  Handle numerical differentiation in MINOS, returning empty g
% 020409  hkh  Adding Mode to Prob.Mode, nState to Prob.nState
% 020616  hkh  No call to nlp_g if NumDiff == 6, because no check in nlp_g
% 050616  hkh  Avoid isfield, assume Prob.FUNCS.fg is defined
% 060814  med  FUNCS used for callbacks instead
% 110724  hkh  Test on NumDiff ~= 6, not < 6, to enable parfor options
