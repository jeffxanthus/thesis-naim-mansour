% Print result of TOMLAB /CPLEX run
%
% function cpxPrint(PriLev,Inform,x,f_k,glnodes,lpiter,Prob,slack,v,basis,rc)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2007 by Tomlab Optimization Inc., $Release: 11.0.0$
% Written Sept 22, 2002.   Last modified Feb 23, 2007.

function cpxPrint(PriLev,Inform,x,f_k,glnodes,lpiter,Prob,slack,v,basis,rc)

global MAX_x MAX_c

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nTOMLAB /CPLEX solved Problem %d:\n',Prob.P);
   fprintf(    '-->-->-->-->-->-->-->-->-->-->\n\n');

   fprintf('Inform = %d. ',Inform);
   [ExitText,ExitFlag] = cplexStatus(Inform);
   fprintf(ExitText);
   fprintf('\n');   
   fprintf('\nObjective function at x (obj) %25.16f --- ',f_k);
   if glnodes > lpiter
      fprintf('Nodes visited%7d. ',glnodes);
   else
      fprintf('LP iterations%7d. ',lpiter);
   end
   fprintf('\n');
   %fprintf('Number of infeasibilities%7d. ',ninf);
   %fprintf('Sum of infeasibilities %e. ',sinf);
   %fprintf('\n');
end

if PriLev > 1
   if isempty(MAX_x)
      MAX_x=length(x);
   end
   fprintf('Optimal x = \n');
   xprinte(x(1:min(length(x),MAX_x)),'x:  ');
end

if PriLev > 2
   fprintf('Slack variables s =\n');
   xprint(slack,'s:');
end

if PriLev > 3
   if isempty(MAX_c)
      MAX_c=20;
   end
   fprintf('Dual variables (Lagrangian multipliers) v = \n');
   xprinte(v(1:min(length(v),MAX_c)),'v:');

   fprintf('Reduced costs r =\n');
   xprint(rc(1:min(length(rc),MAX_c+MAX_x)),'r: ',' %14.9f',5);
end
if PriLev > 4
   fprintf('Basis b =\n');
   xprint(basis(1:min(length(basis),MAX_x)),'b: ',' %14.9f',5);
end

% MODIFICATION LOG:
%
% 020922 hkh  Written
% 021003 ango Cosmetic changes to output
% 070223 hkhh Use cplexStatus to get ExitText from Inform value
