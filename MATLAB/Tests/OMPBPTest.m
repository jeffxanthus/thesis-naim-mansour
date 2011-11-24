function [SNRMatrix SNRorigMatrix] = OMPBPTest(meth,fileName)
%REGULARIZATIONTEST 
addpath('../Samples')
addpath('../')
global methodChoice;
global regularization;
regularization=0;

if nargin < 2
    fileName='BachHymn.wav';
end
method=meth;
methodChoices=[1 3]; 
amountOfClippingLevels=8;
amountOfSamples=4;
clip=[0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2];
length=3000;
ReconstructedMatrix=zeros(3,amountOfClippingLevels,amountOfSamples,length);
% MissingMatrix=zeros(amountOfFiles,amountOfSamples,length);

SNRMatrix=zeros(3,amountOfClippingLevels);
SNRorigMatrix=zeros(3,amountOfClippingLevels);

%Change the input signal, so no false clipping is applied!!!!!!!!!!!!
%Exclude outliers (calculate variance fraction -> compare to value)
[data,largeData,mediumData,smallData,tinyData,fs,noBits] = InitializeTestVariables(fileName,10000);
for k=methodChoices
    methodChoice=k;
    for i=1:amountOfClippingLevels
        input=Clip(tinyData,clip(1,i));   
%         hold off
%         plot(mediumData)
%         hold on
%         plot(input,'r.')
%         pause
        for j=1:amountOfSamples
            disp(['Now simulating for method ' num2str(k) ', clipping level ' num2str(clip(1,i)) ', and sample ' num2str(j)])
            [ReconstructedMatrix(k,i,j,:) dummy]=CSMain(input,method,fs);
            SNR=Evaluation(tinyData,squeeze(ReconstructedMatrix(k,i,j,:))',fs,noBits);
            SNRorig=Evaluation(tinyData,input,fs,noBits);
            SNRMatrix(k,i)=SNRMatrix(k,i)+SNR
            SNRorigMatrix(k,i)=SNRorigMatrix(k,i)+SNRorig
    %         0)),ReconstructedMatrix(i,j,:),fs,noBits);
        end
%         SNRMatrix(t,i)=SNRMatrix(t,i)./length(regParams);
        SNRMatrix(k,i)=SNRMatrix(k,i)./amountOfSamples;
        SNRorigMatrix(k,i)=SNRorigMatrix(k,i)./amountOfSamples;
    end
end

% SNRMatrix=SNRMatrix';
SNRMatrix
SNRorigMatrix

subplot(2,1,1); plot(SNRMatrix(1,:),'.');
subplot(2,1,2); plot(SNRMatrix(2,:),'.');
end

