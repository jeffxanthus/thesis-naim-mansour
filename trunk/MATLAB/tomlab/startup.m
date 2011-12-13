% startup.m
%
% Find and set the paths needed for TOMLAB runs
%
% The user must go through this file and change if 0/if 1 to get
% the correct paths depending on the particular configuration of the system

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2011 by Tomlab Optimization Inc., $Release: 7.7.0$
% Written June 16, 1999.  Last modified May 5, 2011.

echo off

checkarch

tomlabdir = which('tomlablic');
if(isempty(tomlabdir))
   disp(' ');
   disp('---------------------------------');
   disp('ERROR! No Tomlab directory found.');
   disp('---------------------------------');
   disp(' ');
   return
end
TOM=fileparts(tomlabdir);

disp(['The TOMLAB v7.8 directory is ' TOM]);

if isempty(findstr('tomlab',lower(TOM)))
   disp(' ');
   disp('---------------------------------');
   disp('ERROR! No Tomlab directory found.');
   disp('---------------------------------');
   disp(' ');
   disp('Make sure your Tomlab home directory is [...]/tomlab');
   disp(' ');
   return
end

curdir = pwd;
cd(TOM);

BASE=fullfile(cd,'base');
disp(['Found the base path as  ' BASE]);

CGO=fullfile(cd,'cgo');
disp(['Found the cgo path as  ' CGO]);

TP=fullfile(cd,'testprob');
disp(['Found the testprob path as  ' TP]);

EX=fullfile(cd,'examples');
disp(['Found the examples path as  ',EX]);

OPTIM=fullfile(cd,'optim');
disp(['Found the optim path as  ',OPTIM]);

% Change to "if 0" if having the Optimization TB, avoiding TOMLAB versions
% DEFAULT if 1
if 1
   addpath(OPTIM);
   disp(' ');
   disp('Warning - optimization toolbox drop-in replacements');
   disp('Remove tomlab\optim from path to disable');
   disp(' ');
end

MEXINTF=fullfile(cd,'mex');
disp(['Found the mex path as  ',MEXINTF]);

SPLINE=fullfile(cd,'splines');
disp(['Found the splines path as  ',SPLINE]);

QG=fullfile(cd,'quickguide');
disp(['Found the quickguide path as  ',QG]);

MAD=fullfile(cd,'mad');
disp(['Found the MAD path as  ',MAD]);

LIB=fullfile(cd,'lib');
disp(['Found the lib path as  ',LIB]);

USERSG=fullfile(cd,'usersguide');
disp(['Found the usersguide path as  ',USERSG]);

AMPL=fullfile(cd,'ampl');
disp(['Found the ampl path as  ',AMPL]);

COMMON=fullfile(cd,'common');
disp(['Found the common path as  ',COMMON]);

FINANCE=fullfile(cd,'finance');
disp(['Found the finance path as  ',FINANCE]);

addpath(BASE,LIB,CGO,MEXINTF,TP,EX,QG,USERSG,AMPL,COMMON,FINANCE);

MODELLIB=fullfile(cd,'modellib');
disp(['Found the modellib path as  ',MODELLIB]);
MODELLIBPATHS = genpath(MODELLIB);
addpath(MODELLIBPATHS);

% Change to "if 0" if having the Spline TB
%if 1
%%DEFAULTif 1
%   addpath(SPLINE);
%end

% Check if Spline TB routine csapi is in PATH
if ~exist('csapi','file')
   addpath(SPLINE);
end

% Initialize MAD. If you have MAD installed elsewhere, you need to manually
% put the path to the MAD directory in the variable MAD.
% Also, for MATLAB R14 and R2006 users, see mad/@extras/README_EXTRAS.txt

V=version;
if (str2num(V(1:3)) < 6.49)
   fprintf('Warning: Matlab 6.5 or higher needed for MAD\n');
else
   currPath = pwd;
   addpath(MAD);
   cd(MAD)
   startupmad
   cd(currPath);
   clear currPath
end

% Add paths to the main Tomlab directory
addpath(TOM);

% Check if TOMLAB /CPLEX is installed together
% with this installation of TOMLAB
CPLEX=fullfile(cd,'cplex');
if exist(CPLEX,'dir')==7
    addpath(CPLEX);
    addpath(fullfile(CPLEX,'examples'));
    addpath(fullfile(CPLEX,'network'));
end

% Check if TOMLAB /GUROBI is installed
GUROBI=fullfile(cd,'gurobi');
if exist(GUROBI,'dir')==7
    addpath(GUROBI)
    addpath(fullfile(GUROBI,'examples'));
end

% Check if TOMLAB /Xpress is installed together with
% this installation of TOMLAB
XPRESS=fullfile(cd,'xpress');
if exist(XPRESS,'dir')==7
    addpath(XPRESS);
    addpath(fullfile(XPRESS,'examples'));
end

% PROPT
PROPT=fullfile(cd,'propt');
if exist(PROPT,'dir')==7
   disp( ['Found the PROPT path as ' PROPT] );
   PROPTPATHS = genpath(PROPT);
   addpath(PROPTPATHS);
end
clear PROPT PROPTPATHS

% TOMSYM
TOMSYM=fullfile(cd,'tomsym');
if exist(TOMSYM,'dir')==7
   disp( ['Found the TOMSYM path as ' TOMSYM] );
   TOMSYMPATHS = genpath(TOMSYM);
   addpath(TOMSYMPATHS);
end
clear TOMSYM TOMSYMPATHS

MODELLIB=fullfile(cd,'modellib');
disp(['Found the modellib path as  ',MODELLIB]);

% Check if the tomlab\shared directory is set in the Windows PATH
% (LD_LIBRARY_PATH for Unix)

if isunix
   % Unix is case-sensitive. Don't use lower();
   if(strcmp(computer,'MACI'))
     WV='DYLD_LIBRARY_PATH';
   else
     WV='LD_LIBRARY_PATH';
   end
   WP=getenv(WV);
   SH=fullfile(TOM,'shared');
   if isempty(findstr(WP,fullfile('tomlab','shared')))
      % Path to tomlab/shared not found. When setenv exists, try to add it:
      try
         fprintf('\nWarning: the TOMLAB shared files directory %s could not be found in the %s variable.\n',SH,WV);
         fprintf('\nTOMLAB will add ''%s'' to %s using the setenv command.\n',SH,WV);
         setenv(WV,[WP,pathsep,SH]);
      catch
         % setenv missing if really old Matlab, or other problem.
         fprintf('\nWarning: the TOMLAB shared files directory %s could not be found in the %s variable.\n',SH,WV);
         fprintf('Some solvers like CPLEX and Xpress may not work properly. This is purely informational and does not apply\n');
         fprintf('to all operating systems, and also depends on the TOMLAB components that are installed and used.\n');
         fprintf('NOTE: This is not related to the MATLAB path setting, but the OS environment variable %s\n\n',WV);
      end
  end
else
   WP=getenv('PATH');
   SH=fullfile(TOM,'shared');
   if isempty(findstr(lower(WP),fullfile('tomlab','shared')))
      % Path to tomlab/shared not found. When setenv exists, try to add it:
      try
         fprintf('\nWarning: the TOMLAB shared files directory %s could not be found in the Windows PATH variable.\n',SH);
         fprintf('\nTOMLAB will add ''%s'' to PATH using the setenv command.\n',SH);
         setenv('PATH',[WP,pathsep,SH]);
      catch
         % setenv missing if really old Matlab, or other problem.
         fprintf('\nWarning: the TOMLAB shared files directory %s could not be found in the Windows PATH variable.\n',SH);
         fprintf('Some solvers like CPLEX and Xpress may not work properly. This is purely informational and does not apply\n');
         fprintf('to all operating systems, and also depends on the TOMLAB components that are installed and used.\n');
         fprintf('NOTE: This is not related to the MATLAB path setting, but the Windows PATH environment variable.\n\n');
      end
   end
end
clear WP WV SH V v

clear MEXINTF MEXINTF7 TP EX LDO LDODEMO LIB USERSG
clear OPT1X OPTIM SPLINE QG
clear CPLEX GUROBI XPRESS AMPL XA MAD BASE CGO COMMON
clear FINANCE MODELLIB MODELLIBPATHS

s = which('tomlablic','-all');
if length(s) > 2 % There should be 2 (.m and .dll) but not more:
   fprintf(['\nPossible TOMLAB installation error - more than one tomlablic.m is visible.\n' ...
         'This usually means that two TOMLAB installation are visible to Matlab.\n' ...
         'The command ''rehash toolboxcache'' in Matlab might help clear the paths.\n' ...
         'Also, check File... Set Path... in Matlab to make sure that only one\n' ...
         'TOMLAB version is in the Matlab path.\n']);
   clear s curdir tomlabdir
   error('Possible TOMLAB installation error');
else
   clear s
end
disp(' ')
disp('Run tomlablic to display license information');
disp('Done with setting TOMLAB paths');
cd(curdir)
clear tomlabdir

cd(curdir)
clear curdir TOM
