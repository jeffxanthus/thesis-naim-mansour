% pen.m 
%
% Interface to TOMLAB PENSDP and PENBMI solvers
%
% function [x, fx, uoutput, iresults, fresults, info] = pen(p)
%
% --------------------------------------------------------------------
%
% Inputs: 
%
% p        Structure in PENSDP or PENBMI format. If p.ki_dim is present,
%          the problem is assumed to be a bmi problem and is passed to 
%          the MEX file solver penbmi. Otherwise, pensdp is called.
%
% --------------------------------------------------------------------
%
% Outputs: 
%
% x        Optimal point, if found. 
% fx       Function value at x
% uoutput  Vector of multipliers at solution
% iresults Vector of integer results
% fresults Vector of float results
% info     Error flag

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Feb 12, 2003.  Last modified Jan 11, 2005.

function [x, fx, uoutput, iresults, fresults, info] = pen(p, feas)

% feas is still here for compatibility reasons.

if nargin < 2,feas = [];end
if isempty(feas),feas=0;end

% Check for BMI problem
if( isfield(p,'ki_dim') )
  
  [x,info,fx,uoutput,iresults,fresults] = penbmi(p, 1);
  
else
  
   [x,info,fx,uoutput,iresults,fresults] = pensdp(p, 1);

end

% MODIFICATION LOG
% 041214 frhe Changed the calls for the new PENOPT mexes.
% 050113 frhe Readded feas in the argument list.