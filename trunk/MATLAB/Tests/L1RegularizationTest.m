function [SNRMatrix] = L1RegularizationTest(clip)
%REGULARIZATIONTEST 
addpath('../Samples')
addpath('../')
global methodChoice;
global regularization;

method=1;
methodChoice=3;
amountOfFiles=4;
amountOfSamples=4;
regParams=[0.1 0.06 0.04 0.01 0.001]% 0.04 0.02 0.01 0.001]
length=10000;
ReconstructedMatrix=zeros(amountOfFiles,amountOfSamples,length);
% MissingMatrix=zeros(amountOfFiles,amountOfSamples,length);
SNRMatrix=zeros(amountOfFiles,amountOfSamples);
SNRmMatrix=zeros(amountOfFiles,amountOfSamples);

FileNameVector=['BachHymn.wav'; 'BeethoP5.wav';'Forsake2.wav';'Secret_2.wav'];

t=1;
for k=regParams
    regularization=k;
    for i=1:amountOfFiles
        [data,largeData,mediumData,smallData,tinyData,fs,noBits] = InitializeTestVariables(FileNameVector(i,:),10000);
        input=Clip(data,clip);
        for j=1:amountOfSamples
            disp(['Now simulating for parameter ' num2str(k) ', file ' num2str(i) ', and sample ' num2str(j)])
            [ReconstructedMatrix(i,j,:) dummy]=CSMain(input(1,101+(j*20000):101+length-1+(j*20000)),method,fs);
            SNR=Evaluation(data(1,101+(j*20000):101+length-1+(j*20000)),squeeze(ReconstructedMatrix(i,j,:))',fs,noBits);
            SNRMatrix(t,i)=SNRMatrix(t,i)+SNR;
    %         SNRmMatrix(i,j)=Evaluation(data(1,101+(j*20000):10100+(j*2000
    %         0)),ReconstructedMatrix(i,j,:),fs,noBits);
        end
%         SNRMatrix(t,i)=SNRMatrix(t,i)./length(regParams);
        SNRMatrix(t,i)=SNRMatrix(t,i)./5;
    end
    t=t+1;
end

SNRMatrix=SNRMatrix';

subplot(amountOfFiles,1,1); plot(SNRMatrix(1,:),'.');
subplot(amountOfFiles,1,2); plot(SNRMatrix(2,:),'.');
end

