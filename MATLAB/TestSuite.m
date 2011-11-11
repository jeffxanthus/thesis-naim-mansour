function [] = TestSuite()
%TESTSUITE 

dataO=[];
SNRres=[];
SNRmres=[];
k=1;

for i=[0.75 0.85];% 0.75 0.7 0.65 0.6 0.55 0.5];
    i
    disp(['Iteration ' int2str(k) ', for clipping ratio ' num2str(i)]);
    [dataOrig, reconstruction, dummy1, dummy2, SNR, SNRm]=Simulation(i,2,1,0);
    pause
    sound(dataOrig,44100)
    pause
    sound(reconstruction,44100)


    SNRres=[SNRres SNR];
    SNRmres=[SNRmres SNRm];
    k=k+1;
end
figure(2)
subplot(2,1,1);plot(SNRres);
axis([i(1) i(end) 0 40])
subplot(2,1,2);plot(SNRmres);
axis([i(1) i(end) 0 40])
end