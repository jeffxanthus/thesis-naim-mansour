% probAssign implements the TOMLAB (TQ) format for most types of
% optimization.
%
% However, there are a number of specialized routines for different types:
%
% When the type is lp or qp, lpAssign and qpAssign are recommended
% For mip mipAssign must be used.
% For glc glcAssign should be used if there any integer constraints
% For glb either probAssign or glcAssign could be used.
% For ls/cls clsAssign is strongly recommended
% For uc/con conAssign is recommended. The pattern of the Hessian (or
% the constraint gradient matrix) could then be given, and the routine names
% are directly given, without the extra call to tomFiles.
%
% probAssign is setting the variables normally needed for an optimization in
% the TOMLAB structure Prob.
% An additional call to tomFiles.m is needed to set names of functions used:
%    Prob = tomFiles(Prob, .....)
%
% optType is one of the optimization problem types defined in TOMLAB
%
% This routine is usable if you want to call a solver directly and do not
% want to use the TOMLAB Init File format
%
% function Prob = probAssign(optType, x_L, x_U, Name, x_0, fLowBnd, ...
%                 A, b_L, b_U, c_L, c_U, x_min, x_max, f_opt, x_opt);
%
% INPUT (Call with at least four parameters)
%
% optType Any of lp,uc,qp,con,ls,cls,glb,glc
%         For lp,  better to use lpAssign instead
%         For qp,  better to use qpAssign instead
%         For mip, you must use mipAssign instead
%         For glc, you should use glcAssign if there are integer constraints
%         For glb, you could use glcAssign instead
%         For ls/cls,  better to use clsAssign instead
%         For uc/con,  conAssign has expanded functionality
% x_L     Lower bounds on parameters x. If [] set as a nx1 -Inf vector.
% x_U     Upper bounds on parameters x. If [] set as a nx1  Inf vector.
% Name    The name of the problem (string
% x_0     Starting values, default nx1 zero vector
%
% Note:   The number n of the unknown variables x are taken as
%         max(length(x_L),length(x_U),length(x_0))
%         You must specify at least one of these with correct length,
%         then the others are given default values
%
%         The following parameters are optional, and problem type dependent
%         Set empty to get default value
%
% fLowBnd A lower bound on the function value at optimum. Default -1E300
%         A good estimate is not critical. Use [] if not known at all.
%
% A       Matrix A in linear constraints b_L<=A*x<=b_U. Stored dense or sparse.
% b_L     Lower bound vector in linear constraints b_L<=A*x<=b_U.
% b_U     Upper bound vector in linear constraints b_L<=A*x<=b_U.
% c_L     Lower bound vector in nonlinear constraints c_L<=c(x)<=c_U.
% c_U     Upper bound vector in nonlinear constraints c_L<=c(x)<=c_U.
%         c(x) must be defined in a function, see e.g. con_c.m
% x_min   Lower bounds on each x-variable, used for plotting
% x_max   Upper bounds on each x-variable, used for plotting
% f_opt   Optimal function value(s), if known (Stationary points)
% x_opt   The x-values corresponding to the given f_opt, if known.
%         If only one f_opt, give x_opt as a 1 by n vector
%         If several f_opt values, give x_opt as a length(f_opt) by n matrix
%         If adding one extra column n+1 in x_opt, 0 is min, 1 saddle, 2 is max.
%         x_opt and f_opt is used in printouts and plots.
%
%
% Set the variable as empty if this variable is not needed for the particular
% kind of problem you are solving

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.7.0$
% Written Oct 1, 1998.    Last modified Dec 13, 2006.

function Prob = probAssign(optType, x_L, x_U, Name, x_0, fLowBnd, ...
                A, b_L, b_U, c_L, c_U, x_min, x_max, f_opt, x_opt)

if nargin < 15
   x_opt=[];
   if nargin < 14
      f_opt=[];
      if nargin < 13
         x_max=[];
         if nargin < 12
            x_min=[];
            if nargin < 11
               c_U=[];
               if nargin < 10
                  c_L=[];
                  if nargin < 9
                     b_U=[];
                     if nargin < 8
                        b_L=[];
                        if nargin < 7
                           A=[];
                           if nargin < 6
                              fLowBnd=[];
                              if nargin < 5
                                 x_0 = [];
end, end, end, end, end, end, end, end, end, end, end

optType=lower(optType);

n = max([length(x_L),length(x_U),length(x_0)]);

Prob       = ProbDef(1);
Prob = checkAssign(Prob, n, x_0, x_L, x_U, b_L, b_U, A);
Prob.P     = 1;
Prob.N     = n;
Prob.Name  = deblank(Name);
if isempty(x_min)
   Prob.x_min = -1*ones(n,1);
else
   Prob.x_min = x_min;
end
if isempty(x_max)
   Prob.x_max = ones(n,1);
else
   Prob.x_max = x_max;
end
Prob.f_opt = f_opt;
Prob.x_opt = x_opt;

if strcmp(optType,'uc')
   probType = checkType('uc');
elseif strcmp(optType,'qp')
   probType = checkType('qp');
elseif strcmp(optType,'con')
   probType = checkType('con');
elseif strcmp(optType,'ls')
   probType = checkType('ls');
   Prob.LS.SepAlg=0;
   Prob.LS.weightType=0;
   Prob.JacPattern=[];
elseif strcmp(optType,'lls')
   probType = checkType('lls');
   Prob.LS.SepAlg=0;
   Prob.LS.weightType=0;
   Prob.JacPattern=[];
elseif strcmp(optType,'exp')
   probType = checkType('exp');
   Prob.LS.SepAlg=0;
   Prob.LS.weightType=0;
   Prob.JacPattern=[];
elseif strcmp(optType,'cls')
   probType = checkType('cls');
   Prob.LS.SepAlg=0;
   Prob.LS.weightType=0;
   Prob.JacPattern=[];
elseif strcmp(optType,'mip')
   fprintf('To make a MIP problem, use mipAssign\n');
   error('Illegal type of optimzation problem');
elseif strcmp(optType,'lp')
   probType = checkType('lp');
elseif strcmp(optType,'glb')
   probType = checkType('glb');
elseif strcmp(optType,'glc')
   probType = checkType('glc');
   Prob.MIP.IntVars = [];
else
   probType=0;
end

Prob.probType=probType;
Prob.probFile=0;

if ~isempty(fLowBnd) 
   Prob.f_Low=max(Prob.f_Low,fLowBnd); 
end

mN = max(length(c_L),length(c_U));
Prob.mNonLin = mN;

if mN > 0
   if isempty(c_L) 
      Prob.c_L=-Inf*ones(mN,1);
   else
      Prob.c_L=full(double(c_L(:)));
   end
   if isempty(c_U) 
      Prob.c_U=Inf*ones(mN,1);
   else
      Prob.c_U=full(double(c_U(:)));
   end
   if any(Prob.c_L>Prob.c_U)
      error('c_L and c_U have crossover values');
   end
end

% Set Print Level to 0 as default
Prob.PriLevOpt=0;

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print
if isempty(MAX_x)
   MAX_x=20;
end
if isempty(MAX_c)
   MAX_c=20;
end
if isempty(MAX_r)
   MAX_r=30;
end

%Prob.PartSep.pSepFunc= 1;
%Prob.PartSep.index   = 0;

% MODIFICATION LOG
%
% 981012  hkh  f_Low was set wrongly as Prob.f_Low
% 981110  hkh  Change to MAX_x,MAX_c,MAX_r
% 981221  hkh  Must init special global params: Prob.GLOBAL=glbDef;
% 990116  hkh  Set Prob.NLLS.SepAlg=0 if ls, cls or exp type.
% 990117  hkh  Set Prob.NLLS.weightType=0 if ls, cls or exp type.
% 990527  hkh  define glc type
% 990527  mbk  Call glcDef and condef for 'glc' problems.
% 990616  hkh  Call lpDef with two parameters
% 990630  hkh  Add creation of empty mideva.m file.
% 000925  hkh  Change to fLowBnd in LineParam, remove default optParams
% 040102  hkh  Add definition of fields mLin and mNonLin
% 040607  med  Help updates
% 041201  hkh  Add check on number of columns in A, should be n
% 041222  med  Checking lengths for x_L, x_U and x_0
% 050117  med  mlint revision
% 060206  med  Added checks for crossover bounds
% 060705  med  Updated help
% 060818  hkh  Set Prob.f_Low=max(Prob.f_Low,fLowBnd); if fLowBnd~=[]
% 060822  med  All vectors set to full and double
% 061213  med  Moved most input checks to checkAssign