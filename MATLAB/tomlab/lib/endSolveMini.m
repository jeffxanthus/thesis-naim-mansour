% Routine to be called after optimization for LP, QP, LLS etc.
%
% function Result=endSolveMini(Prob,Result)
%
% Computing CPU time and real time elapsed

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written June 7, 2008.   Last modified Aug 13, 2009.

function Result=endSolveMini(Prob,Result)

if isempty(Result)
   return
end

Result.CPUtime = cputime-Prob.TIME0;
Result.REALtime = etime(clock,Prob.TIME1);

% Put Prob into Result struct ALWAYS
Result.Prob=Prob;

% Assign f_k, g_k, r_k, J_k
if ~isempty(Result.x_k)
   x_k = Result.x_k;
   if Result.FuncEv > 0 & isempty(Result.f_k)
      Result.f_k = nlp_f(x_k, Prob);
   end
   if Result.ResEv > 0 & isempty(Result.r_k)
      Result.r_k = nlp_r(x_k, Prob);
   end
   if Result.GradEv > 0 & isempty(Result.g_k)
      Result.g_k = nlp_g(x_k, Prob);
   end
   if Result.JacEv > 0 & isempty(Result.J_k)
      Result.J_k = nlp_J(x_k, Prob);
   end
end

% Add constant f term Prob.fConstant to f_k
Result.f_k = Result.f_k + Prob.fConstant;

% MODIFICATION LOG
%
% 080607 hkh  Written based on endSolve, avoid setting Result.H_k
