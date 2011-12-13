% ====================================================================
% function convflag = CGOisClose(fGoal,f,fTol,nFunc,Iter,PriLev)
% ====================================================================


% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc., $Release: 6.1.0$
% Written April 10, 2008.    Last modified April 10, 2008.

function convflag = CGOisClose(fGoal,f,fTol,nFunc,Iter,PriLev)

convflag = 0;
if isempty(fGoal), return, end
if isinf(fGoal),   return, end

if f <= fGoal
   convflag = 1;
elseif fGoal == 0
   %if abs(f-fGoal) < fTol
   if abs(f) < fTol
      convflag = 2;
   end
elseif abs(f-fGoal) <= abs(fGoal) * fTol
   convflag = 3;
end

if convflag > 0 & PriLev > 0 
   if convflag == 1
      fprintf('Function value %f is less than fGoal %f. ',f,fGoal);
   elseif convflag == 2
      fprintf('Error in function value %f is ',f);
      fprintf('%f <= fTol %f. ',abs(f-fGoal),fTol);
   elseif convflag == 3
      fprintf('Relative error in function value %f is ',f);
      fprintf('%f <= fTol %f. ',abs(f-fGoal)/abs(fGoal),fTol);
   end
   fprintf('# of f(x) evals %d ',nFunc);
   fprintf('and iterations %d\n',Iter);
end


% MODIFICATION LOG:
%
% 080410 hkh Written 
