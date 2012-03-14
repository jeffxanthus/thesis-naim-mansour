function [dataToUse,rec,origSamples, missingSamples,SNR,SNRm,SNROrig,ODG,ODGOrig] = LPFilter(input, reconstruction, ratio)

dataToUse=wavread('BachPartita44k.wav')';
noBits=16;
fs=44000;

[SNROrig ODGOrig]=Evaluation(dataToUse,input,fs); 

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
missingSamples=reconstruction(:,origSamples);
reconstructionOrig=reconstruction;

x=dct(reconstructionOrig, length(reconstructionOrig));
cutOff=round(ratio*length(reconstructionOrig));
xNew=x(1,1:cutOff);
reconstructionOrig=idct(xNew,length(reconstructionOrig));

% LP only on clipped samples
% reconstruction(:,origSamples)=reconstructionOrig(:,origSamples);
reconstruction=reconstructionOrig;

% sound(input,44100)
% pause
% sound(reconstruction,44100)
% pause

subplot(3,1,1);plot(originalSamples,'r.')
title('Clipped values in original signal')
subplot(3,1,2);plot(missingSamples,'.')
title('"Clipped values" in reconstructed signal')
subplot(3,1,3);plot(abs((originalSamples-missingSamples)./originalSamples))
title('Relative error')
axis([0 length(originalSamples) 0 0.5])
% pause
[SNR ODG]=Evaluation(dataToUse,reconstruction,fs,noBits)
SNRm=Evaluation(originalSamples,missingSamples,fs,noBits)
rec=reconstruction;

ODG
pause
end

