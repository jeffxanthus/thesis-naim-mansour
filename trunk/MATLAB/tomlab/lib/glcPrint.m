% Print results during glcFast run
 
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Mar 24, 2002.   Last modified May 25, 2003.

function glcPrint(Iter,nFunc,funMin,fMinEQ,IntVars,xBest,Prob)

if isempty(funMin), funMin = Inf; end

if Iter > 0
   fprintf('Iter:%6d Ev:%6d funMin:%20.10f', Iter, nFunc, funMin);
else
   fprintf('         OK:%6d         f(x):%20.10f', nFunc,funMin);
end

if ~isempty(fMinEQ)
   fprintf(' Eq:%12.10f',fMinEQ);
end

if IntVars > 0
   fprintf(' IV:')
   for i = 1:IntVars
       fprintf(' %d',xBest(i)); 
   end
end

fprintf(' ');
fprintf(' %f',xBest(IntVars+1:end)); 
fprintf('\n');

% MODIFICATION LOG:
%
% 020324  hkh  Written with 3 output variables
% 020326  hkh  Added output of constraint infeasibility and xBest
% 030117  hkh  Change name from fMIn to funMin
% 030525  hkh  Also print all trial points, and if feasible