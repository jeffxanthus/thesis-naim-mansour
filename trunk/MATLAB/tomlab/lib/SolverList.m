% SolverList returns a list of all solvers for a particular solvType
%
% function [SolvList,SolvTypeList] =  SolverList(solvType,LargeScale, Silent)
%
% INPUT:
% SolvType   The TOMLAB solvType number, currently 1-15, or the optType name
%    solvType number  optType name  Description
%            1             uc       Unconstrained Optimization (UC), bounds.
%            2             qp       Quadratic Programming (QP)
%            3             con      Constrained Nonlinear Programming (NLP)
%            4             ls       Nonlinear Least Squares (NLLS), bounds.
%            5             lls      Linear Least Squares (LS)
%            6             cls      Constrained Nonlinear Least Squares
%            7             mip      Mixed-Integer Programming
%            8             lp       Linear Programming
%            9             glb      Global optimization (GO), bounds.
%           10             glc      Global optimization (GO), constraints.
%           11             miqp     Mixed-Integer Quadratic Programming (MIQP)
%           12             minlp    Mixed-Integer Nonlinear Programming (MINLP)
%           13             sdp      Semidefinite Programming (SDP), LMI
%           14             bmi      Linear SDP with BMI constraints
%           15             cgo      Costly Global optimization (CGO)
%
% See checkType for the numbers for the following types
%                          ode      Parameter estimation in ODE (ODEFit)
%                          exp      Exponential sum fitting (ExpFit)
%                          oc       Optimal control (OC)
%                          miqq     Mixed-Integer QP with Quadratic constraints
%                          gp       Geometric Programming Problems
%                          mco      Multi-Criteria Optimization
%                          oc       Optimal control
%                          lcp      Standard Linear Complementarity Problem(LCP)
%                          mcp      Polyhedrally constrained variational
%                                   inequality Problem or
%                                   Mixed Complementarity Problem(MCP)
%                          nts      Nonlinear Time Series
%
% LargeScale 1=Large Scale problems, 0 = Small or medium sized
%
% Silent  If Silent == 1, totally silent
%         If Silent == 0 (default), an information text is displayed
%
% OUTPUT:
% SolvList     A list of all solvers for a particular solvType
% SolvTypeList The solvType numbers corresponding to the elements in SolvList

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Nov 29, 1998.   Last modified Aug 13, 2009.

% Not used now:
%           19             nts      Nonlinear Time Series

function [SolvList,SolvTypeList] =  SolverList(solvType,LargeScale, Silent)

if nargin < 3
   Silent = 0;
   if nargin < 2
      LargeScale = 0;    
      if nargin < 1
         solvType = [];
      end
   end
end

[TomV,os,TV]=tomlabVersion;

SolvList=str2mat(...
     'ucSolve','fmins','fminu' ...
    ,'qpe','qplm','qpSolve','qpBiggs','qpopt','qp','quadprog','qp-minos'...
    ,'qld' ...
    ,'nlpSolve','conSolve','sTrustr','constr','minos','npsol','PDCO'...
    ,'snopt','fmincon'...
    ,'leastsq','lsqnonlin' ...
    ,'lsei', 'lssol' ...
    ,'clsSolve','nlssol' ...
    ,'mipSolve','cutplane','Xpress-MP' ...
    ,'lpSimplex','akarmark','lp','linprog','lpsimp2','lpopt','lp-minos' ...
    ,'glbSolve','ego','glbFast','glcSolve','glcFast','glcCluster' ...
    ,'rbfSolve','bqpd','miqpBB','minlpBB','PENSDP','filterSQP','CPLEX' ...
    ,'PENBMI','PDSCO','nlpsolv' ...
    ,'knitro','conopt','oqnlp','xa','nlpqlp','nlpjob','dfnlp','lgo'...
    ,'DualSolve','slsSolve','sqopt','modfit','socs','dido','coplgp','path'...
    ,'milpSolve','arbfmip');

nToUse = 19;

ODE  = checkType('ode');
EXP  = checkType('exp');
MIQQ = checkType('miqq');
GP   = checkType('gp');
OC   = checkType('oc');
MCP  = checkType('mcp');
MCO  = checkType('mco');

%NTS  = checkType('nts');
%LCP  = checkType('lcp');


% 70 are added to the solvers of less quality
SolvTypeList=[...
      1 1 1    ...
      72 72 2 72 2 2 2 2 ...
      2 ...
      3 3 3 3 3 3 3 ...
      3 3 ...
      4 4 ...
      5 5 ...
      6 6 ...
      7 7 11 ...
      8 78 8 8 78 8 8 ...
      9 15 9 10 10 10 ...
      15 2 11 12 13 3 11 ...
      14 3  73 ...
      3 3 12 11 3,MCO, 6 10 ...
      8 6 2,ODE,OC,OC,GP,MCP ...
      7 15];

whichTV=[1 1 1 1 1 1 1 2 1 1 2 1 1 1 1 1 2 4 1 4 1 1 1 1 3 1 3 1 1 8 ...
         1 1 1 1 1 2 2 ...
         1 5 1 1 1 1 ...
         5 7 7 7 6 7 9 ...
         10 1 1 ...
         11 12 14 15 16 16 16 17 ...
         1 1 4 16 24 19 25 20 ...
         1 5];

% How to update with new optType
%
% Add new 2-3 letter code to optType
% Add full name to optTypeDescr
% Add Shorter description to optTypeDescrShort
% Set nToUse to the number of types to actually display
% 
% New solvers are added to SolveList 

optType=str2mat('uc','qp','con','ls','lls','cls', 'mip','lp','glb', ...
       'glc','miqp','minlp','sdp','bmi','cgo','','','','','','','','','','');

optType(ODE,:)  = 'ode  ';
optType(EXP,:)  = 'exp  ';
optType(MIQQ,:) = 'miqq ';
optType(MCO,:)  = 'mco  ';
optType(GP,:)   = 'gp   ';
optType(MCP,:)  = 'mcp  ';
optType(OC,:)   = 'oc   ';

optTypeDescr=str2mat('Unconstrained Optimization (UC), bounds', ...
    'Quadratic Programming (QP)', 'Constrained Nonlinear Programming (NLP)', ...
    'Nonlinear Least Squares (NLLS), bounds', 'Linear Least Squares (LS)', ...
    'Constrained Nonlinear Least Squares', 'Mixed-Integer Programming', ...
    'Linear Programming', 'Global optimization (GO), bounds', ...
    'Global optimization (GO), constraints', ...
    'Mixed-Integer Quadratic Programming (MIQP)', ...
    'Mixed-Integer Nonlinear Programming (MINLP)', ...
    'Semidefinite Programming (SDP)', 'Linear SDP with BMI constraints', ...
    'Costly Global optimization (CGO)', ...
    '', '', '', '','','','','');

optTypeDescr(ODE,:)  = 'Parameter estimation in ODE (ODEFit)       ';
optTypeDescr(EXP,:)  = 'Exponential sum fitting (ExpFit)           ';
optTypeDescr(MIQQ,:) = 'Mixed-Integer QP w Quadratic constr (MIQQ) ';
optTypeDescr(MCO,:)  = 'Multi-Criteria Optimization (MCO)          ';
optTypeDescr(GP,:)   = 'Geometric Programming (GP)                 ';
optTypeDescr(MCP,:)  = 'Var InEq or Mixed Complementarity   (mcp)  ';
optTypeDescr(OC,:)   = 'Optimal control (OC)                       ';

optTypeDescrShort=str2mat('UC', 'QP', 'NLP', 'NLLS', 'LS', ...
    'Constrained Nonlinear NLLS', 'Mixed-Integer Programming', ...
    'Linear Programming', 'GO, bounds', 'GO, constraints', ...
    'MIQP', 'MINLP', 'SDP /LMI', 'Linear SDP /BMI', ...
    'CGO', '','','', '','','','','');

optTypeDescrShort(ODE,:)  = 'ODEFit                    ';
optTypeDescrShort(EXP,:)  = 'ExpFit                    ';
optTypeDescrShort(MIQQ,:) = 'MIQQ                      ';
optTypeDescrShort(MCO,:)  = 'Multi-Criteria Opt        ';
optTypeDescrShort(GP,:)   = 'Geometric Programming     ';
optTypeDescrShort(MCP,:)  = 'LCP, MCP                  ';
optTypeDescrShort(OC,:)   = 'Optimal Control           ';


if isempty(solvType) % No input argument was given
%   SolvList = SolvList(SolvTypeList < 70,:);
   if ~Silent

      SolvList = [];
      % Print information of possible input arguments
      fprintf('\nPlease give the TOMLAB optType name ');
      fprintf('or the solvType number, currently 1-9, as input\n\n');
      fprintf('optType  solvType  Description\n');
      for t = 1:nToUse
          fprintf(optType(t,:));
         for s = 1:(9-size(optType(t,:),2))
            fprintf(' ');
         end
         if t < 10
            fprintf('%d         ',t);
         else
            fprintf('%d        ',t);   
         end
         fprintf(optTypeDescr(t,:));
         fprintf('\n');
      end
      fprintf('\n');
   end
   return;
end

i = [];
if ischar(solvType)
   % Find recommended solver
   RecSolvLarge=GetSolver(solvType,1);
   RecSolvSmall=GetSolver(solvType,0);
   % Convert string to the solvType number
   i=strmatch(deblank(lower(solvType)),optType,'exact');
   solvType=i;
   if isempty(solvType)
      % Given solvType was given as a string but it did not match
      % any of 'uc','qp','con','ls','exp','cls', 'mip','lp','glb','glc',
      % See checkType for correct list
      SolvList = [];
      SolvTypeList = [];
      return;
   end
else
   solvType = min(nToUse,max(1,solvType));
   solvString = deblank(optType(solvType(1),:));
   %RecSolvLarge=[];
   %RecSolvSmall=[];
   %for i=1:size(solvstring,1)
   %    RecSolvLarge=[RecSolveLarge,GetSolver(solvString(i,:),1)];
   %    RecSolvSmall=[RecSolveSmall,GetSolver(solvString(i,:),0)];
   %end
   RecSolvLarge=GetSolver(solvString,1);
   RecSolvSmall=GetSolver(solvString,0);
end

% Read probTypeList that gives the list of what other types of problem
% a solver of type solvType will handle

ix = [];
for i=1:nToUse
    if i ~= solvType
       [DataFile,NameFile,DefFile,probTypeList]=nameprob(i,0);
       CanSolve = any(solvType==probTypeList);
       if CanSolve
          % Add the type i to the list of the types that can solve solvType
          ix = [ix i];
       end
   end
end

% Make a list of solvers of the types in ix,
% also able to solve solvType problem
iopt = [];
for i = 1:size(ix,2)
   iopt=cat(2,iopt,find(ix(i)==SolvTypeList));
end
OthSolvList=SolvList(iopt,:);
OthSolvTV=whichTV(iopt);

% Make a list of the solvers for the particular solvType
ipart=find(solvType==SolvTypeList);
PartSolvList=SolvList(ipart,:);
PartSolvTypeList=SolvTypeList(ipart);
PartSolvTV=whichTV(ipart);

% Check which ones there is a license for
% Check TV vector

SolvList = PartSolvList;
SolvTypeList = PartSolvTypeList;
if Silent, return; end

% Print recommended choices of solvers

if LargeScale
   fprintf('\nTomlab recommended choice for large scale ');
   fprintf(optTypeDescr(solvType,:));
   fprintf('\n\n');
   fprintf(RecSolvLarge);
   fprintf('\n');
end
fprintf('\nTomlab recommended choice for ');
if LargeScale
   fprintf('small scale ');
end
fprintf(optTypeDescr(solvType,:));
fprintf('\n\n');
fprintf(RecSolvSmall);
fprintf('\n\n');

fprintf('Other solvers for ');
fprintf(deblank(optTypeDescrShort(solvType,:)));
fprintf('\n\n');

fprintf('   Licensed:\n\n   ');
count = 0;
for i = 1:size(PartSolvList)   
   if TV(PartSolvTV(i)) & ...
   ~(strcmpi(deblank(PartSolvList(i,:)),RecSolvSmall) | ...
    (strcmpi(deblank(PartSolvList(i,:)),RecSolvLarge) &...
      LargeScale))
      fprintf(PartSolvList(i,:));
      fprintf('\n   ');
      count = 1;
   end
end
if ~count
   fprintf('NONE\n');
end
fprintf('\n');

fprintf('   Non-licensed:\n\n   ');
count = 0;
for i = 1:size(PartSolvList)   
   if ~TV(PartSolvTV(i))& ...
   ~(strcmpi(deblank(PartSolvList(i,:)),RecSolvSmall) | ...
    (strcmpi(deblank(PartSolvList(i,:)),RecSolvLarge) &...
      LargeScale))
      fprintf(PartSolvList(i,:));
      fprintf('\n   ');
      count = 1;
   end
end
if ~count
   fprintf('NONE\n');
end
fprintf('\n');

% Print other solvers also capable of solving problem type

fprintf('Solvers also handling ');
fprintf(optTypeDescrShort(solvType,:));
fprintf('\n\n');

fprintf('   Licensed:\n\n   ');
count = 0;
for i = 1:size(OthSolvList)   
   if TV(OthSolvTV(i))
      fprintf(OthSolvList(i,:));
      fprintf('\n   ');
      count = 1;
   end
end
if ~count
   fprintf('NONE\n');
end
fprintf('\n');

fprintf('   Non-licensed:\n\n   ');
count = 0;
for i = 1:size(OthSolvList)   
   if ~TV(OthSolvTV(i))
      fprintf(OthSolvList(i,:));
      fprintf('\n   ');
      count = 1;
   end
end
if ~count
   fprintf('NONE\n');
end
fprintf('\n');

% MODIFICATION LOG
%
% 981129  hkh  Written
% 981203  mbk  Check if isempty(solvType) then solvType was given as input
%              argument but it did not match any of 'uc','qp', ...
% 981209  hkh  Add lpsimp2 among LP solvers
% 990906  hkh  Add glcSolve and glcRun
% 990913  hkh  Add MIP as number 7, and use use exp, not ef for number 5.
% 001004  hkh  Adding lpopt
% 010714  hkh  Adding xpress-mp and qp-xpress-mp
% 010715  hkh  Adding glbFast
% 010815  hkh  Adding glcFast
% 011111  hkh  Adding glcCluster and rbfSolve
% 020701  hkh  Adding six new solvers
% 020702  hkh  Update optType with new problem types
% 020708  bjo  Changed output to recommended solvers + other solvers.
%              Also prints which are licensed and which are not.
% 030117  hkh  Change type miqq to bmi
% 030123  hkh  Correct the list whichTV, add PDCO, PDSCO
% 030213  ango Correct names (Dundee) and filterSQP type
% 040126  hkh  Now only returns the solvers for the particular type
% 040126  hkh  Added flag Silent, if true totally silent
% 040126  hkh  Empty argument displays all available solvers
% 040413  med  Added more solvers
% 040414  hkh  Correct ego, now type glc (10)
% 050502  hkh  Add CGO,ODE,OC. Only display nToUse many
% 050602  hkh  Add GP. Make more general solution, using checkType
% 050801  med  isstr replaced by ischar
% 060212  hkh  Add arbfmip
% 090813  med  mlint check