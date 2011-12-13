% optim_dc.m
%
% function dc = optim_dc(x, Prob, varargin)
%
% optim_dc is used to implement the OPT TB 2.0 interface
%
% The constraint gradient dc is returned,
% if available in the global variable NLP_dc
% with the corresponding x value in NLP_xdc
%
% optim_dc is called from the TOMLAB gateway function nlp_dc.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written July 29, 1999.  Last modified Jan 14, 2003.

function dc = optim_dc(x, Prob)

global NLP_xdc NLP_dc

% Return the constraint gradient, if computed.

if ~isempty(NLP_xdc)
   if all(x==NLP_xdc)
      dc=NLP_dc;
   elseif ~isempty(NLP_dc)
      c = optim_cdc(x, Prob);
      dc=NLP_dc;
   else
      dc=[];
   end
else
   dc=[];
end

% MODIFICATION LOG:
%
% 001218 hkh Return one row for each constraint, tomlab format
% 030115 hkh Minor revision