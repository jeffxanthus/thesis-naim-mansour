% Initialization of structure Prob
%
% function Prob = ProbDef(New)
%
% INPUT:
%
% New     Flag if a new type of problem is defined.
%         If New ~=0, initialize probType as empty
%         Also set counters:
%              n_f n_g n_H n_c n_dc n_d2c n_J n_r n_d2r
%         as double(0)
%
% OUTPUT:
%  Prob   Structure

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written May 18, 1998.   Last modified July 23, 2011.

function Prob = ProbDef(New)

global probType n_f n_g n_H n_c n_dc n_d2c n_J n_r n_d2r 

if nargin < 1
   New=[];
end

if isempty(New), New=0; end

if New % Init counters to 0
   probType = [];
   n_f=double(0); n_g=double(0); n_H=double(0); n_c=double(0); n_dc=double(0); 
   n_d2c=double(0); n_r=double(0); n_J=double(0); n_d2r=double(0);
end

LP = LineParamDef;
FUNCS = struct('f',[], 'g',[], 'H',[], 'c',[], 'dc',[], 'd2c',[], 'r',[], ...
              'J',[], 'd2r',[],'fc',[],'gdc',[],'fg',[],'cdc',[], 'rJ',[]);

SOL = struct( 'SpecsFile', [], 'PrintFile', [], 'SummFile', [], ...
    'xs', [], 'hs', [], 'nS', double(0), 'hElastic', [], ...
    'iState', [], 'cLamda', [], 'R', [], 'optPar', -999*ones(1,71), ...
    'optParN', 71);

DUNDEE = struct('optPar', -999*ones(1,20));
PENOPT = struct('ioptions', -999*ones(1,12), 'foptions', -999*ones(1,12));
Solver = struct('Alg',0,'Name',[],'Method',0);

% Find out Matlab version
v = version;
y = str2num(v(1:1));
z = str2num(v(3:4));

Prob = struct('TOMLAB','v7.8','MATLAB',[y,z], 'Threads', double(0), 'A', [], ...
  'ADObj', double(0), 'ADCons', double(0), 'BIG', [], 'b_L', [], 'b_U', [],...
  'c_L', [], 'c_U', [], 'CheckNaN', 0, 'cName', [], 'cols', [], 'ConIx', [],...
  'ConsDiff',double(0), 'ConsIdx', [], 'ConsPattern', [], 'd2cPattern', [],...
  'd2LPattern',[],'FAST', 0, 'fConstant',double(0),'f_Low',-1E300, ...
  'GATEF',0,'GATEC',0, ...
  'f_opt', [], 'g_k',[],'GradTolg', [], 'GradTolH', [], 'GradTolJ', [],...
  'HessIx',[],'HessPattern', [], 'JacIx', [],'JacPattern', [], ...
  'LargeScale', [], 'MaxCPU', inf ,'MENU', 0, ...
  'Mode', 2, 'nState',double(1), 'N', [], 'mLin',[],'mNonLin',[], ...
  'Name','User Problem 1', 'NumDiff', double(0), 'P', double(1),...
  'plotLine', 0, 'PriLev', double(0), 'PriLevOpt', double(0),  ...
  'probFile', [], 'probType', [], 'rows',[], 'simType', double(0), ...
  'smallA', double(0), ...
  'SolverDLP', [], 'SolverFP', [], 'SolverLP', [], 'SolverQP', [], ...
  'uP', [], 'uPName', [], 'WarmStart',0,'Warning',1, 'x_0', [], 'x_L', [],...
  'x_U', [], 'x_min', [], 'x_max', [], 'x_opt', [], 'xName', [], ...
  'QP', [], 'LS', [], 'MIP', [], 'GO', [], 'CGO', [], 'ExpFit', [], ...
  'NTS', [], 'LineParam', LP, 'optParam', [], 'PartSep', [], ...
  'Solver',Solver,'FUNCS',FUNCS, 'CONOPT',[], 'DUNDEE',DUNDEE, 'GENO',[], ...
  'glbDirect',[], 'glcDirect',[], 'GP',[], 'KNITRO',[], 'LGO',[], ...
  'LSGRG2',[], 'MILPSOLVE',[], 'MISQP',[], 'OQNLP',[], 'PENOPT',PENOPT, ...
  'SOL', SOL, 'CHECK',0);

% MODIFICATION LOG:
%
% 980825  hkh  Changed nts to NTS for Nonlinear Time Series.
% 980909  mbk  Structure field GLOBAL added.
% 980920  hkh  Change Prob.f_min to Prob.optParam.f_Low. Delete Prob.optPar
% 980921  mbk  optparam.f_Low changed to optParam.f_Low.
% 981006  hkh  Added name uPName, defined together with uP to avoid conflicts
% 981010  hkh  Changed function call logic, using TOMLAB gateway routines
% 981018  hkh  New field NLLS. t and Yt moved to NLLS field.
% 981020  hkh  New field PartSep, for partially separable functions
% 981023  hkh  Added Prob.p_d2r and Prob.FUNCS.p_d2r
% 981026  hkh  Add field NumDiff for numerical differentiation
%              Changed back f_Low to top level of struct
% 981105  hkh  Delete field f_0
% 981111  hkh  Add definition of Solver field in Prob
% 981119  hkh  Add input variable New.
% 990213  hkh  prob.x_min=[]; prob.x_max=[]; should be capital P
% 000710  hkh  Add SOL field parameters
% 000928  hkh  Remove GLOBAL field
% 001014  hkh  Added field ConsDiff
% 001019  hkh  Always use fixed names for the TOMLAB gateway routines
% 001105  hkh  Make field NLLS into LS, add linear least squares
% 010903  hkh  Now 63 parameters in optPar
% 011110  hkh  Add two new fields, GO and RBF, used for global optimization
% 020103  hkh  Change field RBF to CGO
% 020409  hkh  Add field Mode = request for functions and gradients
% 020409  hkh  Add field nState,= 1 if 1st time call for functions & gradients
% 020512  hkh  Add field DUNDEE for Fletcher and Leyffer Dundee solvers
% 020630  hkh  Add field PENSDP for Tomlab /PENSDP
% 030107  hkh  Change field PENSDP to PENOPT, for /PENSDP and /PENBMI
% 030117  hkh  Add field CheckNaN, if NaN in derivatives should be checked
% 030129  hkh  Add field Tomlab with version number
% 030522  hkh  Add field fc and gdc in FUNCS
% 030524  hkh  Add field simType, default 0, Use double(1), double(0)
% 031129  hkh  Add fields ADObj, ADCons, use double for nState and counters
% 040102  hkh  Add fields mLin mNonLin for lengths of constraints
% 040106  hkh  Add fields ConIx, ConsIdx, and rearrange in alphabetic order
% 040106  hkh  Add fields JacIx, HessIx
% 040412  hkh  Add field Prob.Warning, put x_L and x_U after each other
% 040425  hkh  Add field Prob.smallA and Prob.MaxCPU
% 040608  hkh  Improve comments
% 041213  hkh  Add field BIG, default []
% 050422  hkh  Add field fConstant, default 0
% 050605  hkh  Change to optParN = 65, for SNOPT 6
% 050606  hkh  Change to optParN = 71, for SNOPT 7
% 050616  hkh  Add fields FUNCS.fg and FUNCS.cdc, avoid isfield in /SOL callbacks
% 050616  hkh  Add field g_k, avoid isfield in nlp_g
% 050725  med  DUNDEE and PENOPT added optPar (options) vector default
% 050902  med  Added d2cPattern (used by KNITRO and CONOPT)
% 051216  med  Add field FUNCS.rJ, avoid isfield in NLSSOL calls
% 060814  med  FUNCS used for callbacks instead
% 080311  hkh  NumDiff, ConsDiff set [], not 0
% 080416  hkh  Add [] fields for CONOPT, GENO, glbDirect, glcDirect, GP;
% 080416  hkh  and KNITRO, LGO, LSGRG2, MILPSOLVE, OQNLP. Also MISQP.
% 080604  hkh  Try to set NumDiff=ConsDiff=0
% 080606  med  Minor clean ups
% 080608  hkh  Add field FAST; default 0. GATEF, GATEC used in gateways
% 080919  med  Removed extra LineParamDef
% 081117  med  Version to 7.0 with TomSym platform
% 091001  hkh  Add field MATLAB with Matlab version in numeric form
% 100910  ango MATLAB field now 2 elements, [major,minor] due to 7.10 and higher
% 101124  ango TOMLAB version.
% 110723  hkh  Field Threads for Parallel Computing Toolbox, using parfor
