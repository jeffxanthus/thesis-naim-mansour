% checkType returns true if the input type (probType) is the same
% as the wished type given as a string (strType)
%
% It may also return the string name of the type of optimization, 
% given the corresponding number, or the number given the string name of 
% the type of optimization 
%
% function  isType = checkType(strType,numType)   (1) 
% OR       strType = checkType(numType)           (2) 
% OR       numType = checkType(strType)           (3) 
%
% INPUT: 
%
% strType    The problem type given as a string
% numType    The problem type given as a number
%
% OUTPUT: 
%
% Mode 1: isType = checkType(strType,numType) 
%
% isType  True if problem type defined by strType is the same as by numType
%         If isempty(strType), the string Type corresponding to numType(1) is 
%         given
%         If isempty(numType), the number corresponding to strType is given
%
% Example: 
%         isType = checkType('con',1);  gives isType == 0  (false)
%         isType = checkType('con',3);  gives isType == 1  (true)
%
% Mode 2: strType = checkType(numType)
%
% strType The character representation of the given numerical problem
%         type.
% 
% Example: 
%         strType = checkType(3);   gives strType = 'con';
%
% Mode 3: numType = checkType(strType)
%
% Example:
%         numType = checkType('con'); gives numType = 3
%
%
% Current types of optimization problems in TOMLAB:
%
%      1.  uc    Unconstrained Optimization
%      2.  qp    Quadratic Programming
%      3.  con   Constrained Optimization (nonlinear programming)
%      4.  ls    Nonlinear Least Squares
%      5.  lls   Linear Least Squares
%      6.  cls   Constrained Nonlinear Least Squares
%      7.  mip   Mixed-Integer Programming
%      8.  lp    Linear Programming
%      9.  glb   Box-bounded Global Optimization
%      10. glc   Constrained Global Optimization, also integer variables
%      11. miqp  Mixed-Integer Quadratic Programming (MIQP)
%      12. minlp Mixed-Integer Nonlinear Programming (MINLP)
%      13. sdp   Semidefinite Programming, Linear SDP with LMI constraints
%      14. bmi   Linear SDP with BMI Constraints
%      15. cgo   Costly Global Optimization (constraints & integer variables)
%      16. ode   Parameter estimation in ODEs
%      17. exp   Parameter estimation in exponential models
%      18. miqq  Mixed-Integer Quadratic Programming with Quadratic 
%                constraints
%      19. gp    Geometric Programming Problems
%      20. mco   Multi-Criteria Optimization
%      21. oc    Optimal control
%      22. lcp   Standard Linear Complementarity Problem (LCP)
%      23. mcp   Polyhedrally constrained variational inequality Problem or
%                Mixed Complementarity Problem(MCP)
%      24. nts   Nonlinear Time Series

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Nov 5, 1998.    Last modified Jun 2, 2005.

function isType = checkType(varargin)

if nargin == 2
   strType = varargin{1};
   numType = varargin{2};
elseif nargin == 1
   if ischar(varargin{1})
      strType = varargin{1};
      numType = [];
   else
      strType = [];
      numType = varargin{1};
   end
elseif nargin < 1
   error('checkType needs at least one input argument')
end

if isempty(numType)
   if nargin == 2
      isType = 0;
      return
   end
   switch lower(strType)
      case 'uc' 
         isType=1;
      case 'qp' 
         isType=2;
      case 'con' 
         isType=3;
      case 'ls' 
         isType=4;
      case 'lls' 
         isType=5;
      case 'cls' 
         isType=6;
      case 'mip' 
         isType=7;
      case 'lp' 
         isType=8;
      case 'glb' 
         isType=9;
      case 'glc' 
         isType=10;
      case 'miqp' 
         isType=11;
      case 'minlp' 
         isType=12;
      case 'sdp' 
         isType=13;
      case 'bmi' 
         isType=14;
      case 'cgo' 
         isType=15;
      case 'ode' 
         isType=16;
      case 'exp' 
         isType=17;
      case 'miqq' 	
         isType=18;   
      case 'gp' 	
         isType=19;
      case 'mco' 	
         isType=20;   
      case 'oc' 
         isType=21;
      case 'lcp' 	
         isType=22;
      case 'mcp' 	
         isType=23;
      case 'nts' 
         isType=24;
      otherwise
         isType=[];
   end
   return
end

for i = 1:length(numType)
   switch numType(i)
      case 1 
         Type='uc';
      case 2
         Type='qp';
      case 3
         Type='con';
      case 4
         Type='ls';
      case 5
         Type='lls';
      case 6
         Type='cls';
      case 7
         Type='mip';
      case 8
         Type='lp';
      case 9 
         Type='glb';
      case 10
         Type='glc';
      case 11
         Type='miqp';
      case 12
         Type='minlp';
      case 13
         Type='sdp';
      case 14
         Type='bmi';
      case 15
         Type='cgo';
      case 16
         Type='ode';
      case 17
         Type='exp';
      case 18 	
         Type='miqq';   
      case 19 	
         Type='gp';   
      case 20
         Type='mco';   
      case 21
         Type='oc';
      case 22 	
         Type='lcp';
      case 23 	
         Type='mcp';
      case 24
         Type='nts';
      otherwise
         Type=' ';
   end

   if isempty(strType)
      isType = Type;
      return
   end

   isType = strcmpi(Type,strType);
 
   if isType, return; end
end

% MODIFICATION LOG
%
% 001105 hkh Written
% 020701 hkh Adding four new types, miqp, minlp, sdp, miqq.
% 021010 hkh Set isType=0 if probType is empty (may occur in tomGUI)
% 030117 hkh Change miqpp to bmi
% 040419 med Added lcp and mcp
% 040517 med Added miqq
% 050503 hkh Added ode, cgo, nts, moved some (exp, miqq)
% 050601 med gp added (24)
% 050601 ang Allows one arg (string or number) or two args (string & number)
% 050602 hkh Switch numbers for gp and others
