function checkdll(full)
% checkdll calls the DLLs (or libs) to check the license
%
% function checkdll(full)
%
% If full >= 1 check all solver MEX
%
% Default full == 1

if nargin < 1
    full = 1;
end

disp('.............................................................')
disp('-------------------------------------------------------------')
disp('checkdll calls the MEX files to check the installation')
disp('--------------------------')
disp('Test of license file for Tomlab')
disp('--------------------------')

try
   tomlablic
catch
   fprintf('Error while checking tomlablic.\n\n');
   GetEnvironment;

   try
      tomlablic(80)
   catch
      fprintf('Failed to execute tomlic/tomlablic(80)')
   end
   % Rethrow last error
   rethrow(lasterror);
end

[TomV,os,TV] = tomlabVersion;
if full & TV(1)
    pause(2)
    disp('--------------------------')
    disp('Test of Tfzero solver license')
    disp('--------------------------')
    checkMex('Tfzero');
    pause(2)
    disp('------------------------------------------')
    disp('Test of Tomlab lsqr (Tlsqr) solver license')
    disp('------------------------------------------')
    checkMex('Tlsqr');
    pause(2)
    disp('--------------------------')
    disp('Test of LSEI solver license')
    disp('--------------------------')
    checkMex('lsei');
    pause(2)
    disp('--------------------------')
    disp('Test of TNNLS solver license')
    disp('--------------------------')
    checkMex('Tnnls');
    pause(2)
    disp('--------------------------')
    disp('Test of TOMSOL solver license')
    disp('--------------------------')
    checkMex('tomsol');
end
if TV(2)
    disp('--------------------------')
    disp('Test of MINOS solver license')
    disp('--------------------------')
    checkMex('minos');
elseif TV(1)
    disp('--------------------------')
    disp('Test of QLD solver license')
    disp('--------------------------')
    checkMex('qld');
end
if full & TV(2)
    pause(2)
    disp('--------------------------')
    disp('Test of QPOPT solver license')
    disp('--------------------------')
    checkMex('qpopt');
    pause(2)
    disp('--------------------------')
    disp('Test of LPOPT solver license')
    disp('--------------------------')
    checkMex('lpopt');
end
if full & TV(4)
    pause(2)
    disp('--------------------------')
    disp('Test of SQOPT solver license')
    disp('--------------------------')
    checkMex('sqopt');
    pause(2)
    disp('--------------------------')
    disp('Test of SNOPT solver license')
    disp('--------------------------')
    checkMex('snopt');
end
if full & TV(3)
    pause(2)
    disp('--------------------------')
    disp('Test of NPSOL solver license')
    disp('--------------------------')
    checkMex('npsol');
    pause(2)
    disp('--------------------------')
    disp('Test of LSSOL solver license')
    disp('--------------------------')
    checkMex('lssol');
    disp('--------------------------')
    disp('Test of NLSSOL solver license')
    disp('--------------------------')
    checkMex('nlssol');
end
if full & TV(5)
    pause(2)
    disp('--------------------------')
    disp('Test of CGO Toolbox solver license')
    disp('--------------------------')
    checkMex('tomsol',27);
end
if full & TV(6)
    pause(2)
    disp('--------------------------')
    disp('Test of PENSDP solver license')
    disp('--------------------------')
    checkMex('pensdp');
end
if full & TV(7)
    pause(2)
    disp('--------------------------')
    disp('Test of MINLP solver licenses, sparse and dense versions')
    disp('--------------------------')
    checkMex('bqpds');
    checkMex('bqpdd');
    checkMex('miqpBBs');
    checkMex('miqpBBd');
    checkMex('minlpBBs');
    checkMex('minlpBBd');
    checkMex('filSQPs');
    checkMex('filSQPd');
end
if full & TV(8)
    pause(2)
    disp('--------------------------')
    disp('Test of Xpress interface license');
    disp('--------------------------')
    checkMex('xpressmp');
end
if full & TV(9)
    pause(2)
    disp('--------------------------')
    disp('Test of CPLEX interface license');
    disp('--------------------------')
    checkMex('cplexmex');
end
if full & TV(10)
    pause(2)
    disp('--------------------------')
    disp('Test of PENBMI solver license');
    disp('--------------------------')
    checkMex('penbmi');
    checkMex('penbmiQ');
end
if full & TV(11)
    pause(2)
    disp('--------------------------')
    disp('Test of KNITRO solver license');
    disp('--------------------------')
    checkMex('Tknitro');
end
if full & TV(12)
    pause(2)
    disp('--------------------------')
    disp('Test of CONOPT solver license');
    disp('--------------------------')
    checkMex('conopt');
end
if full & TV(13)
    pause(2)
    disp('--------------------------')
    disp('Test of AMPL interface licenses');
    disp('--------------------------')
    checkMex('amplqp');
    checkMex('amplfunc');
    checkMex('spamfunc');
end
if full & TV(14)
    pause(2)
    disp('--------------------------')
    disp('Test of OQNLP interface licenses');
    disp('--------------------------')
    checkMex('oqnlp');
end
if full & TV(22)
    pause(2)
    disp('--------------------------')
    disp('Test of MSNLP interface licenses');
    disp('--------------------------')
    checkMex('msnlp');
    checkMex('lsgrg2');
end

disp('.............................................................')

disp('.............................................................')
disp('If everything works you can put a comment on the last row in')
disp('startup.m (the call to checkdll) to avoid this test')
disp('every time you start Tomlab')


function checkMex(name,varargin)

if exist(name)==3
    try
       feval(name,varargin{:})
    catch
       fprintf('Error while checking %s.\n\n',name);
       GetEnvironment;
       rethrow(lasterror);
       pause
    end
else
    disp(['Warning: Licensed MEX file ' deblank(name) ' is not properly installed. Please check your setup.']);
end


% GetEnvironment displays useful information about the Matlab and OS
% environment for debugging purposes

function GetEnvironment

cc = '';
cv = '';

v = version;

if ispc
   [i,w] = dos('ver');
   w(find(w==13))=' ';
   p = getenv('PATH');
else
   [i,w]     = unix('uname -a');

   % Try to find gcc version (and version)
   [i,cc]    = unix('which gcc');
   if ~isempty(cc)
      [i,cv] = unix([cc ' -v']);
   end
   
   if(strcmp(computer,'MAC'))
      p = getenv('DYLD_LIBRARY_PATH');
   else
      p = getenv('LD_LIBRARY_PATH');
   end
end

w   = deblank(w);
lic = license;
tl  = getenv('TOMLAB_LICENSE_FILE');

fprintf('\n\n=========================================\n');
fprintf('System environment: %s\nML: %s\n\n',deblank(w),lic);
   
if ~isempty(tl)
   fprintf('TOMLAB_LICENSE_FILE=%s\n\n',tl);
else
   fprintf('TOMLAB_LICENSE_FILE is not set\n\n');
end

if ispc
   fprintf('PATH: %s\n',p);
else
   if(strcmp(computer,'MAC'))
      fprintf('DYLD_LIBRARY_PATH: %s\n',p);
   else
      fprintf('LD_LIBRARY_PATH: %s\n',p);
   end
end

fprintf('=========================================\n\n');
