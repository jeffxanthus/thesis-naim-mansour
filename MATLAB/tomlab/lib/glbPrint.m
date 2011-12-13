% Print results during glbFast run

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Apr 11, 2002.   Last modified Jan 7, 2003.

function glbPrint(Iter,nFunc,funMin,xBest,Prob)

if isempty(funMin), funMin = Inf; end

fprintf('Iter:%6d Ev:%6d funMin:%20.10f', Iter, nFunc, funMin);

fprintf(' ');
fprintf(' %f',xBest(IntVars+1:end));
fprintf('\n');

% MODIFICATION LOG:
%
% 020411  hkh  Written with 3 output variables, and xBest
% 030116  hkh  Use name funMin instead of fMin