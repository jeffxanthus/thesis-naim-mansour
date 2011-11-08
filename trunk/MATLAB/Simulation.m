function [dataToUse,reconstruction,SNR,SNRm] = Simulation(clip,size,method,clipamount)
%TESTSUITE Summary of this function goes here
%   clip: clipping percentage (between 0 & 1)
%   size: 1=small, 2=medium, 3=large
%   method: 1=> Ax=y, 2=>Ax+Be=y

[data, largeData, mediumData, smallData, tinyData, fs, noBits]=InitializeTestVariables;

switch size
    case {0}
        dataToUse=tinyData;
    case {1}
        dataToUse=smallData;
    case {2}
        dataToUse=mediumData;
    case {3}
        dataToUse=largeData;
end

if nargin > 3
    input=Clip(dataToUse,clip,clipamount);
else
    input=Clip(dataToUse,clip);
end

[reconstruction, missingSamples]=CSMain(input,method,fs);

origSamples=[];
testData=dataToUse(1:length(reconstruction),1)';

Mn=[];
Mp=[];
input=input(1,1:length(reconstruction));
MaxS=max(input);
MinS=min(input);
for t=1:length(input) 
    if(input(1,t)>=MaxS || input(1,t)<=MinS)...
        && ((t==1 || t==length(input)) || (input(1,t-1)-input(1,t+1))==0 || (((length(find(MaxS)))==1) && input(1,t)==MaxS)...
        || (((length(find(MinS)))==1) && input(1,t)==MinS))
        origSamples=[origSamples t];
    end
end
origSamples=testData(:,origSamples);

subplot(3,1,1);plot(origSamples,'.')
title('Clipped values in original signal')
subplot(3,1,2);plot(missingSamples,'.')
title('"Clipped values" in reconstructed signal')
subplot(3,1,3);plot(abs((origSamples-missingSamples)./origSamples))
title('Relative error')
axis([0 length(origSamples) 0 0.5])
pause
SNR=Evaluation(dataToUse(1:length(reconstruction),1),reconstruction,fs,noBits)
SNRm=Evaluation(5*origSamples',5*missingSamples',fs,noBits)
end

