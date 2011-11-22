function [SNRMatrix] = AxBeTest(mChoice,fileName)
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
amountOfClippingLevels=3;
amountOfSamples=4;
clip=[0.9 0.8 0.7];% 0.6 0.5 0.4 0.3 0.2];
length=10000;
ReconstructedMatrix=zeros(2,amountOfSamples,length);
% MissingMatrix=zeros(amountOfFiles,amountOfSamples,length);
SNRMatrix=zeros(2,amountOfSamples);
SNRmMatrix=zeros(2,amountOfSamples);

[data,largeData,mediumData,smallData,tinyData,fs,noBits] = InitializeTestVariables(fileName,10000);
for k=method
    for i=1:amountOfClippingLevels
        input=Clip(data,clip(1,i));
        for j=1:amountOfSamples
            disp(['Now simulating for method ' num2str(k) ', clipping level ' num2str(clip(1,i)) ', and sample ' num2str(j)])
            [ReconstructedMatrix(i,j,:) dummy]=CSMain(input(1,101+(j*20000):101+length-1+(j*20000)),k,fs);
            SNR=Evaluation(data(1,101+(j*20000):101+length-1+(j*20000)),squeeze(ReconstructedMatrix(i,j,:))',fs,noBits);
            SNRMatrix(k,i)=SNRMatrix(k,i)+SNR;
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

