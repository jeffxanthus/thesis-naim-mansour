function [data,largeData,mediumData,smallData,tinyData,fs,noBits] = InitializeTestVariables()
%INITIALIZETESTVARIABLES 

[data,fs,noBits]=wavread('bach_partita.wav');

data=data(:,1);
largeData=data(10000:2000000,1);
mediumData=data(100000:250000,1);
smallData=data(20000:30000,1);
tinyData=data(10000:11000,1);
end

