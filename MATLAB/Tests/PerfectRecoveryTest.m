function [] = PerfectRecoveryTest(sparsity,length)
%PERFECTRECOVERYTEST Summary of this function goes here
%   Detailed explanation goes here
global methodChoice
global regularization

close all;
methodChoice=3;
regularization=0;
[freq data]=SparseSignalConstructor(sparsity,length);
method=1;
fs=44100;
minimumSamples=2*((2*sparsity-1).^2); %NO NOISE YET -- THEN factor 4
minimumSamples
length
if minimumSamples>length
    disp('No perfect recovery possible')
%     return;
end
k=1;
% minimumSamples=length-100;
sound(data,fs)
for i=1002
    disp(['Iteration ' int2str(k) ', for amount of samples ' int2str(i)]);
    input=Clip(data,950/1000);
    subplot(2,1,1);plot(data)
    subplot(2,1,2);plot(input)
    pause
    MaxI=max(input);
    MinI=min(input);
    size(find(input>=MaxI))
    size(find(input<=MinI))
%     [reconstruction, missingSamples]=CSMain(input,method,fs);
    [spect reconstruction]=CSDeclip(input);
    
    subplot(2,2,1);plot(data);
    axis([0 length min(data)-1 max(data)+1])
    subplot(2,2,2);plot(reconstruction);
    axis([0 length min(data)-1 max(data)+1])
    subplot(2,2,3);plot(freq);
    axis([0 length min(freq)-1 max(freq)+1])
    subplot(2,2,4);plot(spect);
    axis([0 length min(freq)-1 max(freq)+1])
    size(data)
    size(reconstruction)
    SNR=Evaluation(data',reconstruction,fs)
    SNR2=Evaluation(data',input,fs)
%     SNRres=[SNRres SNR];
    sound(reconstruction,fs)
    pause
    k=k+1;
end


end

