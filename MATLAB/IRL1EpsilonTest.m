function [SNRs SNRms] = IRL1EpsilonTest(clip,size,method)
%IRL1-EPSILONTEST 
% Author: Naim Mansour
global fFf
global methodChoice

methodChoice=1;

epsilons=[0.9 0.95] %1 1.1 1.2 1.3];

SNRs=zeros(length(epsilons),3);
SNRms=zeros(length(epsilons),3);
for i=1:length(epsilons)
    fFf=epsilons(i);
    disp(['The value of epsilon is now ' num2str(fFf)])
    ra=round((rand()+1)*100000);
    c=1;
    for j=[ra-10000 ra ra+10000]
        [dataOrig, reconstruction, dummy1, dummy2, SNR, SNRm]=Simulation(clip,size,method,ra);
        SNRs(i,c)=SNR;
        SNRms(i,c)=SNRm;
        c=c+1;
    end
end

subplot(2,1,1);
hold on;
plot(SNRs(:,1),'r')
plot(SNRs(:,2),'b')
plot(SNRs(:,3),'g')
hold off;
subplot(2,1,2);
plot(SNRms(:,1),'r')
hold on;
plot(SNRms(:,2),'b')
plot(SNRms(:,3),'g')
end

