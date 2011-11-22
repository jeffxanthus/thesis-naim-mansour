function [dataToUse,reconstruction,origSamples, missingSamples,SNR,SNRm,SNROrig,SNRmorig] = Simulation(clip,size,method,offSet,clipamount)
%TESTSUITE Summary of this function goes here
%   clip: clipping percentage (between 0 & 1)
%   size: 1=small, 2=medium, 3=large
%   method: 1=> Ax=y, 2=>Ax+Be=y

[data, largeData, mediumData, smallData, tinyData, fs, noBits]=InitializeTestVariables('BeethoP5.wav',offSet);

fs=44000;
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
if nargin > 4
    input=Clip(dataToUse,clip,clipamount);
else
    input=Clip(dataToUse,clip);
end

SNROrig=Evaluation(dataToUse,input,fs); 

[reconstruction, missingSamples]=CSMain(input,method,fs);

origSamples=[];
MaxS=max(input);
MinS=min(input);
for t=1:length(input) 
    if(input(1,t)>=MaxS || input(1,t)<=MinS)...
        && (((t==1 && input(1,t+1)-input(1,t)==0)  || (t==length(input) && input(1,t)-input(1,t-1)==0))...  
        || ((t~=1 && t~=length(input)) && input(1,t-1)-input(1,t+1))==0)
        origSamples=[origSamples t];
    end
end
originalSamples=dataToUse(:,origSamples);
originalClipped=input(:,origSamples);

subplot(3,1,1);plot(originalSamples,'r.')
title('Clipped values in original signal')
subplot(3,1,2);plot(missingSamples,'.')
title('"Clipped values" in reconstructed signal')
subplot(3,1,3);plot(abs((originalSamples-missingSamples)./originalSamples))
title('Relative error')
axis([0 length(originalSamples) 0 0.5])
% pause
SNR=Evaluation(dataToUse,reconstruction,fs,noBits)
SNRmorig=Evaluation(originalSamples,originalClipped,fs,noBits)
SNRm=Evaluation(originalSamples,missingSamples,fs,noBits)
end

