% Setup the paths to run TOMLAB /CPLEX. Assume we are in the cplex
% directory. If cplexmex is not found, issue an error message

CPX = pwd;

if exist('cplexmex')~=3
  error(['cplex/startupcplex.m: No MEX file cplexmex is present. Please' ...
	 ' install it and re-run startupcplex.m in the cplex directory.']);
else
  disp(['The TOMLAB /CPLEX directory is          ' CPX]);
  addpath(CPX);
end

% Examples directory
EXP=fullfile(CPX, 'examples');
if exist(EXP)==7
  disp(['The TOMLAB /CPLEX Examples directory is ' EXP]);
  
  addpath(EXP);
else
  warning('cplex/startupcplex.m: cplex/examples directory missing. Please check your installation.');
end
  
% Network directory
NET=fullfile(CPX, 'network');
if exist(NET)==7
   disp(['The TOMLAB /CPLEX Network directory is  ' NET]);
   addpath(NET);
else
   warning('cplex/startupcplex.m: cplex/network directory missing. Please check your installation.');
end

% Common directory
olddir = pwd;
cd('..');
COMMON=fullfile(pwd, ['common']);
if exist(COMMON)==7
  disp(['The common directory is                 ' COMMON]);
  addpath(COMMON);
else
  warning('cplex/startupcplex.m: common directory missing. Please check your installation.');
end
cd(olddir);

if exist('tomlablic')==3
  % The license file looks to be installed right here, in the
  % cplex folder, or it is already in the Matlab path. 
  % No action necessary. 
  tomlablic
else
  % Try one level higher. 

  olddir = cd;

  cd('..');
  
  if exist('tomlablic')==3
    LICDIR = pwd;
    disp(['The license directory - the Tomlab directory - is  ' LICDIR]);
    addpath(LICDIR);
    tomlablic
  else
    % Cannot find tomlablic - give up and tell the user
    cd(olddir);
    error(['cplex/startupcplex.m: Cannot find tomlablic, neither in the' ...
	   ' cplex directory nor one level higher (the tomlab' ...
	   ' directory). Please install a valid license file and try' ...
	   ' again']);
  end
  cd(olddir);
end

clear CPX EXP LICDIR NET COMMON olddir