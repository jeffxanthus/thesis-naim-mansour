% Routine for checking preset values in Prob and set all undefined values 
% in Prob. 
%
% function Prob = odeProbCheck(Prob);
%
% INPUT:
%  Prob      Problem structure
%  odeSolver odesolver name
%
% OUTPUT:
%  Prob     Problem structure
           
% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 21, 2005.  Last modified Jul 6, 2005.

function Prob = odeProbCheck(Prob)

% Check the struct and define any undefined items
if isfield(Prob,'ODE')         
   if isfield(Prob.ODE,'CHECK')
      if Prob.ODE.CHECK==1
         % odeProbCheck has already checked this structure
         return
      end
   end
   if ~isfield(Prob.ODE,'P'),       Prob.ODE.P=[]; end
   if ~isfield(Prob.ODE,'absTol'),  Prob.ODE.absTol=[]; end
   if ~isfield(Prob.ODE,'relTol'),  Prob.ODE.relTol=[]; end
   if ~isfield(Prob.ODE,'f'),       Prob.ODE.f=[]; end
   if ~isfield(Prob.ODE,'J'),       Prob.ODE.J=[]; end
   if ~isfield(Prob.ODE,'InitStep'),Prob.ODE.InitStep=[]; end
   if ~isfield(Prob.ODE,'Name'),    Prob.ODE.Name=[]; end
   if ~isfield(Prob.ODE,'tInit'),   Prob.ODE.tInit=[]; end
   if ~isfield(Prob.ODE,'tStop'),   Prob.ODE.tStop=[]; end
   if ~isfield(Prob.ODE,'tWant'),   Prob.ODE.tWant=[]; end
   if ~isfield(Prob.ODE,'Y0'),      Prob.ODE.Y0=[]; end
   if ~isfield(Prob.ODE,'LSODE'),   Prob.ODE.LSODE=[]; end
   if ~isfield(Prob.ODE,'RKSUITE'), Prob.ODE.RKSUITE=[]; end  
   if ~isfield(Prob.ODE,'ML'),      Prob.ODE.ML=[]; end  
   Prob.ODE.CHECK=1;
else
   Prob.ODE=struct('P', double(1), 'absTol',[], 'relTol',[], 'f',[], 'J',[],...
                   'InitStep',[],'Name',[],'tInit',[],'tStop',[],'tWant',[],...
                   'Y0',[], 'LSODE',[], 'RKSUITE',[], 'ML',[], 'CHECK',1);
end

if isempty(Prob.ODE.P),      Prob.ODE.P = double(1);end
if isempty(Prob.ODE.relTol), Prob.ODE.relTol=1E-4; end
if isempty(Prob.ODE.absTol), Prob.ODE.absTol=1E-8; end

% MODIFICATION LOG:
%
% 050421  bkh  Created from copy of ProbCheck
% 050422  bkh  Changed odeH_s/hStart, odeT_s/tStart and odeT_e/tEnd to 
%              InitStep, tInit and tStop respectively
% 050502  hkh  Clean up, remove unnecessary code
% 050502  hkh  Make more efficient. Change order of apperance in ODE structure
% 050502  hkh  Add Prob.ODE.CHECK, if set, no further check is done
% 050502  hkh  Add More default fields that need to set; Prob.ODE.P = 1;
% 050503  hkh  Remove odeSolver from input
% 050705  med  Help updated