%function WSInfo = saveCGO(Mode, CGOsolver, Name, O, F, X, F_m, F00, Cc,...
%                          nInit, Fpen, fMinIdx, rngState)
%
% Mode = 1
% Does nothing, returns WSInfo = [];
%
% Mode = 1
% Saves warm start info on file cgoSave1.mat every iteration
% Copies cgoSave1.mat to cgoSave.mat
% Returns WSInfo = [];
%
% Mode = 2
% Save warm start info on file cgoSave.mat after the run
% Delete cgoSave1.mat if it exists
% Creates warm start structure WSInfo as output
%
% 

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written July 14, 2008.   Last modified Oct 30, 2009.

function WSInfo = saveCGO(Mode, CGOsolver, Name, O, F, X, F_m, F00, Cc,...
                          nInit, Fpen, fMinIdx, rngState)
if Mode == 0
   WSInfo = [];
elseif Mode == 1
   WSInfo = [];
   try
     save('cgoSave1.mat','Name','O','F','X','F_m','F00','Cc',...
          'nInit','Fpen','fMinIdx', 'rngState');
     if exist('cgoSave1.mat','file')
        copyfile('cgoSave1.mat','cgoSave.mat','f');
        %if ~isunix
        %   xx=system('copy cgoSave1.mat cgoSave.mat > cgoSave.out');
        %   copyfile('cgoSave1.mat','cgoSave.mat','f');
        %else
        %   copyfile('cgoSave1.mat','cgoSave.mat','f');
        %end
     end
   catch
     warning('Failed to save intermediate warm start info to cgoSave1.mat');
     err = lasterror;
     disp(err.message);
   end
else
   try
     save('cgoSave.mat','Name','O','F','X','F_m','F00','Cc',...
          'nInit','Fpen','fMinIdx', 'rngState');
     if exist('cgoSave1.mat','file') > 0
        %clear cgoSave1.mat
        try
          delete('cgoSave1.mat');
        catch 
          err = lasterror;
          disp(err.message);
          fprintf('%s: Not possible to delete cgoSave1.mat\n',CGOsolver);
        end
     end
     if exist('cgoSave.out','file') > 0
        %clear cgoSave.out
        try
          delete('cgoSave.out');
        catch
          err = lasterror;
          disp(err.message);
          fprintf('%s: Not possible to delete cgoSave.out\n',CGOsolver);
        end
     end
   catch 
     warning('Failed to save warm start information to cgoSave.mat');
     err = lasterror;
     disp(err.message);
     fprintf('Check that the file is writable and ');
     fprintf('that %s is executing in a writable folder',CGOsolver);
     disp('Warm start information is returned in Result.CGO.WarmStartInfo');
   end

   % Create struct output with warmstart information
   WSInfo = struct(...
     'Name',Name,...
     'O',O,...
     'F',F,...
     'X',X,...
     'F_m',F_m,...
     'F00',F00,...
     'Cc',Cc,...
     'nInit',nInit,...
     'Fpen',Fpen,...
     'fMinIdx',fMinIdx,...
     'rngState',rngState);
end

% MODIFICATION LOG
%
% 081118 hkh  Use copyfile instead of system(copy...
% 090824 hkh  Minor mlint revision
% 091020 hkh  Avoid catch ME, use lasterror
