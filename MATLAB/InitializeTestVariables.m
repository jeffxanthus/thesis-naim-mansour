function [data,largeData,mediumData,smallData,tinyData,fs,noBits] = InitializeTestVariables(offSet)
%INITIALIZETESTVARIABLES 
% Author: Naim Mansour
[data,fs,noBits]=wavread('bach_partita.wav');

data=data(:,1);
largeData=data(10000:2000000,1);
mediumData=data(100000:250000,1);
smallData=data(offSet+20000:offSet+30000,1);
tinyData=data(offSet+10000:offSet+11000,1);
end

