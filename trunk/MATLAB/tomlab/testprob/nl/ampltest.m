% ampltest

function ampltest()

files = dir('*.nl');

nrfiles = size(files,1);

for(i = 1:nrfiles)
%  disp(['Testing AMPL: ' files(i).name]);
  sparse = mod(i,2);
  amplAssign(files(i).name, sparse, [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], 0);
end