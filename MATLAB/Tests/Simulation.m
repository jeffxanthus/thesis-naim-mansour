function [input,reconstruction,origSamples, missingSamples,SNR,SNRm,SNROrig,SNRmorig,ODGmorig,ODGm] = Simulation(clip,size,method,offSet,clipamount)
%TESTSUITE Summary of this function goes here
%   clip: clipping percentage (between 0 & 1)
%   size: 1=small, 2=medium, 3=large
%   method: 1=> Ax=y, 2=>Ax+Be=y

[data, largeData, mediumData, smallData, tinyData, fs, noBits]=InitializeTestVariables('bach_partita.wav',offSet);

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

estimateSparsity(dataToUse)

if nargin > 4
    input=Clip(dataToUse,clip,clipamount);
else
    input=Clip(dataToUse,clip);
end

[SNROrig ODGOrig]=Evaluation(dataToUse,input,fs); 

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

subplot(4,1,1);plot(dataToUse,'r.')
title('Clipped values in original signal')
subplot(4,1,2);plot(input,'g.')
title('Clipped values')
subplot(4,1,3);plot(reconstruction,'.')
title('"Clipped values" in reconstructed signal')
subplot(4,1,4);plot(abs((dataToUse-reconstruction)./dataToUse))
title('Relative error')
% axis([0 length(originalSamples) 0 0.5])
% subplot(6,1,5);plot(dataToUse,'r')
% title('Original signal')
% subplot(6,1,6);plot(reconstruction,'b')
% title('Reconstructed signal')
% pause
[SNR ODG]=Evaluation(dataToUse,reconstruction,fs,noBits)
[SNRmorig ODGmorig]=Evaluation(originalSamples,originalClipped,fs,noBits)
[SNRm ODGm]=Evaluation(originalSamples,missingSamples,fs,noBits)

ODG
ODGOrig
end

