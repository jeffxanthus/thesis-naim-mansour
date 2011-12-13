% copyexfiles 
%
% Script for copying files needed by this example.
% 

files  = {'DefPar.m',
          'cplex.m',
          'cpx2cbinfo.m',
          'cpx2retvec.m',
          'cplexmex'};

nrfiles = length(files);

for(i = 1:nrfiles)
    disp(['Copying file: ' files{i}]);
    status = copyfile(which(files{i}), '.');
    if(status == 0)
        disp(['Failed to copy file: ' files{i}]);
    end
end