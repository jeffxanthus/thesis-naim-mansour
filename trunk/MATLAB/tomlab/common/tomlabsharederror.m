% function tomlabsharederror
%
% Executed when installation error due to shared folder not seen
% by operating system

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.1.0$
% Written Sep 22, 2005.   Last modified Feb 7, 2006.

function tomlabsharederror

% Try to generate the actual path to the expected shared folder
try
   W = fullfile(fileparts(which('tomlablic')),'shared');
catch
   W = fullfile('tomlab','shared');
end


if ispc
  P = 'PATH';
else
  c = computer;
  if strcmp(c,'MAC') | strcmp(c,'MACI')
    P = 'DYLD_LIBRARY_PATH';
  else
    P = 'LD_LIBRARY_PATH';
  end
end

    
fprintf('\nWarning - Possible installation problem. \n \n');
fprintf('%s should be included in the environment variable %s.\n',W,P);
fprintf('Please check getenv(''%s'') in Matlab - it should contain\n',P);
fprintf('%s. Please reboot the computer - if still not \n',W);
fprintf('working add/edit the %s environment variable manually. \n \n',P);
if ~ispc
   fprintf('For other operating systems than Windows refer to README.TOMLAB \n');
   fprintf('included with the license files. \n \n');
end

error('tomlab/shared folder not found');

% MODIFICATION LOG
%
% 050922  med   Written
% 060207  ango  Little more information provided to user. 

