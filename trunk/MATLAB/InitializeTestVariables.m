function [data,largeData,mediumData,smallData,tinyData,fs,noBits] = InitializeTestVariables(offSet)
%INITIALIZETESTVARIABLES 
% Author: Naim Mansour
[data,fs,noBits]=wavread('bach_partita.wav');

data=data(:,1)';
largeData=data(1,10000:2000000);
mediumData=data(1,100000:250000);
smallData=data(1,offSet+20000:offSet+30000);
tinyData=data(1,offSet+10000:offSet+11000);
end

