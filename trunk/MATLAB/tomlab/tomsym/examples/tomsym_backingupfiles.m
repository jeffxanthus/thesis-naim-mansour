%% Backing Up Files
%
%% Problem description
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
%% Variables
%
%  maxuse                     A maximal number of disks
%  capacity                   Storage available on each disk
%  sizes                      Filesizes
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
maxuse   = 10;
capacity = 1440;
sizes    = [46;55;62;87;108;114;137;164;253;364;372;388;406;432;461;851];

f = length(sizes);
d = maxuse;

save = tom('save',f,d,'int');
use  = tom('use',d,1,'int');

% All variables are binary
bnds = {0 <= save <= 1, 0 <= use <= 1};

% An item must be stored on one unit
con1 = {sum(save,2) == 1};

% Stay below capacity of storage unit
con2 = {(sizes'*save)' <= capacity*use};

% Objective
objective = sum(use);

constraints = {bnds, con1, con2};
options = struct;
options.solver = 'cplex';
options.name   = 'Backing up files';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    files        = length(sizes);
    filepos      = sol.save;
    disks_in_use = sum(sol.use);
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
% 090308 med   Converted to tomSym