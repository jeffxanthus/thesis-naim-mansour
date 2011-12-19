function [input reconstruction SNRres SNRmres] = TestSuite()
%TESTSUITE 
global methodChoice
global regularization
global clip
global fL

addpath('../');
dataO=[];
SNRres=[];
SNRmres=[];
k=1;

methodChoice=4;
regularization=0;
fL=30;

for i=[0.55];% 0.75   0.7 0.65 0.6 0.55 0.5];
    i
    fL
    clip=i;
    disp(['Iteration ' int2str(k) ', for clipping ratio ' num2str(i)]);
    [input, reconstruction, dummy1, dummy2, SNR, SNRm,SNROrig,SNRmorig,ODGmorig,ODGm]=Simulation(i,1,1,300000);
    SNROrig
    SNR
    SNRmorig
    SNRm
    pause
    sound(input,44100)
    pause
    sound(reconstruction,44100)
    pause

    SNRres=[SNRres SNR];
    SNRmres=[SNRmres SNRm];
    k=k+1;
end
% figure(2)
% subplot(2,1,1);plot(SNRres);
% axis([i(1) i(end) 0 40])
% subplot(2,1,2);plot(SNRmres);
% axis([i(1) i(end) 0 40])
end