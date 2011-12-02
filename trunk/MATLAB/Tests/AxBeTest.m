function [ReconstructedMatrix SNRMatrix] = AxBeTest(mChoice,fileName)
%REGULARIZATIONTEST 
addpath('../Samples')
addpath('../')
global methodChoice;
global regularization;
regularization=0;

if nargin < 2
    fileName='BachHymn.wav';
end
method=[1 2];
methodChoice=mChoice;
amountOfClippingLevels=8;
amountOfSamples=1;
clip=[0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2];
length=3000;
ReconstructedMatrix=zeros(2,amountOfClippingLevels,amountOfSamples,length);
% MissingMatrix=zeros(amountOfFiles,amountOfSamples,length);
SNRMatrix=zeros(2,amountOfClippingLevels);
SNRmMatrix=zeros(2,amountOfClippingLevels);

%Change the input signal, so no false clipping is applied!!!!!!!!!!!!
%Exclude outliers (calculate variance fraction -> compare to value)
[data,largeData,mediumData,smallData,tinyData,fs,noBits] = InitializeTestVariables(fileName,10000);
for k=method
    for i=1:amountOfClippingLevels
        input=Clip(smallData,clip(1,i));       
        for j=1:amountOfSamples
            disp(['Now simulating for method ' num2str(k) ', clipping level ' num2str(clip(1,i)) ', and sample ' num2str(j)])
            [ReconstructedMatrix(k,i,j,:) dummy]=CSMain(input(1,101+(j*2000):101+length-1+(j*2000)),k,fs);
            SNR=Evaluation(smallData(1,101+(j*2000):101+length-1+(j*2000)),squeeze(ReconstructedMatrix(k,i,j,:))',fs,noBits);
            SNRMatrix(k,i)=SNRMatrix(k,i)+SNR
    %         SNRmMatrix(i,j)=Evaluation(data(1,101+(j*20000):10100+(j*2000
    %         0)),ReconstructedMatrix(i,j,:),fs,noBits);
        end
%         SNRMatrix(t,i)=SNRMatrix(t,i)./length(regParams);
        SNRMatrix(k,i)=SNRMatrix(k,i)./amountOfSamples;
    end
end

% SNRMatrix=SNRMatrix';
SNRMatrix

subplot(2,1,1); plot(SNRMatrix(1,:),'.');
subplot(2,1,2); plot(SNRMatrix(2,:),'.');
end

