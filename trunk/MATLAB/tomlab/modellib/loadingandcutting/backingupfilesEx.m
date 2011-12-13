% function Result = backingupfilesEx(PriLev)
%
% Creates a TOMLAB MIP problem for backing up files
%
% BACKING UP FILES
%
% Before leaving on holiday, you wish to backup your most important
% files onto floppy disks. You have got empty disks of 1.44Mb
% capacity. The sixteen files you would like to save have the
% following sizes: 46kb, 55kb, 62kb, 87kb, 108kb, 114kb, 137kb,
% 164kb, 253kb, 364kb, 372kb, 388kb, 406kb, 432kb, 461kb, and 851kb.
% 
% Assuming that you do not have any program at hand to compress the
% files and that you have got a sufficient number of floppy disks to
% save everything, how should the files be distributed in order to
% minimize the number of floppy disks used?      
%
% VARIABLES
%
% maxuse                     A maximal number of disks
% capacity                   Storage available on each disk
% sizes                      Filesizes
%
% RESULTS
%
% Results are explained by running:
% Result       = backingupfilesEx(2); 
%
% REFERENCES
%
% Applications of optimization... Gueret, Prins, Seveaux
% http://web.univ-ubs.fr/lester/~sevaux/pl/index.html
%
% INPUT PARAMETERS
% PriLev       Print Level
%
% OUTPUT PARAMETERS
% Result       Result structure.

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.2.0$
% Written Oct 10, 2005.   Last modified Mar 14, 2006.

function Result = backingupfilesEx(PriLev)

if nargin < 1
   PriLev = 1;
end

maxuse   = 10;
capacity = 1440;
sizes    = [46;55;62;87;108;114;137;164;253;364;372;388;406;432;461;851];

Prob = backingupfiles(maxuse, capacity, sizes);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   files        = length(sizes);
   filepos      = reshape(Result.x_k(1:files*maxuse),files,maxuse);
   disks_in_use = sum(Result.x_k(files*maxuse+1:files*maxuse+maxuse));
   disp(['a total of ' num2str(disks_in_use) ' disks are in use'])
   for d = 1:maxuse,
      files_on_disk = filepos(:,d);         % what files on disk d
      if (abs(sum(files_on_disk)) >= 0.5),         % disk d is not empty?
         file_ids = find(files_on_disk);
         disp(['  files on one of them: '  num2str(file_ids') ])
      end
   end  
end

% MODIFICATION LOG
%
% 051010 med   Created.
% 060112 per   Added documentation
% 060125 per   Moved disp to end
% 060314 med   checking for empty disc changed