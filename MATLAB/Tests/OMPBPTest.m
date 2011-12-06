function [SNRMatrix SNRorigMatrix ODGMatrix ODGorigMatrix] = OMPBPTest(meth,fileName,frameLength)
%OMPBP 
addpath('../Samples')
addpath('../')
global methodChoice;
global regularization;
global clip;
global fL

regularization=0;

if nargin < 2
    fileName='BachHymn.wav';
end
method=meth;
methodChoices=[5]; 
amountOfClippingLevels=8;
amountOfSamples=4;
clipping=[0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2];
length=10000;
% ReconstructedMatrix=zeros(3,amountOfClippingLevels,amountOfSamples,length);
fL=frameLength;

SNRMatrix=zeros(3,amountOfClippingLevels);
SNRorigMatrix=zeros(3,amountOfClippingLevels);
ODGMatrix=zeros(3,amountOfClippingLevels);
ODGorigMatrix=zeros(3,amountOfClippingLevels);


%Change the input signal, so no false clipping is applied!!!!!!!!!!!!
%Exclude outliers (calculate variance fraction -> compare to value)
t=1;
for k=methodChoices
    methodChoice=k;
    for i=1:amountOfClippingLevels
         clip=clipping(1,i);
        for j=1:amountOfSamples
            [data,largeData,mediumData,smallData,tinyData,fs,noBits] = InitializeTestVariables(fileName,200000+j*100000);
            input=Clip(smallData,clipping(1,i)); 
            [rsa csa]=size(input);
            origSamples=[];
            MaxS=max(input);
            MinS=min(input);
            for d=1:csa 
                if(input(1,d)>=MaxS || input(1,d)<=MinS)...
                    && (((d==1 && input(1,d+1)-input(1,d)==0)  || (d==csa && input(1,d)-input(1,d-1)==0))...  
                    || ((d~=1 && d~=csa) && input(1,d-1)-input(1,d+1))==0)
                    origSamples=[origSamples d];
                end
            end
            originalSamples=smallData(:,origSamples);
            originalClipped=input(:,origSamples);
            [rsb csb]=size(input);
            ReconstructedSamples=zeros(1,csb);
            disp(['Now simulating for method ' num2str(k) ', clipping level ' num2str(clipping(1,i)) ', and sample ' num2str(j)])
            [ReconstructedSignal ReconstructedSamples]=CSMain(input,method,fs);
            
            [dummy1 ODG]=Evaluation(smallData,ReconstructedSignal,fs,noBits)
            [dummy2 ODGorig]=Evaluation(smallData,input,fs,noBits)
           
            [SNR dummy3]=Evaluation(originalSamples,ReconstructedSamples,fs,noBits);
            [SNRorig dummy4]=Evaluation(originalSamples,originalClipped,fs,noBits);
            
            SNRMatrix(t,i)=SNRMatrix(t,i)+SNR
            ODGMatrix(t,i)=ODGMatrix(t,i)+ODG
            SNRorigMatrix(t,i)=SNRorigMatrix(t,i)+SNRorig
            ODGorigMatrix(t,i)=ODGorigMatrix(t,i)+ODGorig
    %         0)),ReconstructedMatrix(i,j,:),fs,noBits);
        end
%         SNRMatrix(t,i)=SNRMatrix(t,i)./length(regParams);
        SNRMatrix(t,i)=SNRMatrix(t,i)./amountOfSamples;
        SNRorigMatrix(t,i)=SNRorigMatrix(t,i)./amountOfSamples;
        ODGMatrix(t,i)=ODGMatrix(t,i)./amountOfSamples;
        ODGorigMatrix(t,i)=ODGorigMatrix(t,i)./amountOfSamples;
    end
    t=t+1;
end

% SNRMatrix=SNRMatrix';
SNRMatrix
SNRorigMatrix
ODGMatrix
ODGorigMatrix

subplot(2,1,1); plot(SNRMatrix(1,:),'.');
subplot(2,1,2); plot(SNRMatrix(2,:),'.');
end

