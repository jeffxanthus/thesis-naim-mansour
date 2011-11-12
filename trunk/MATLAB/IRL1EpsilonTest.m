function [SNRs SNRms] = IRL1EpsilonTest(clip,size,method)
%IRL1-EPSILONTEST 
global fFf

epsilons=[0.8 0.9] %1 1.1 1.2 1.3];

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

subplot(2,1,1);plot(SNRs(:,1),'r',SNRs(:,2),'b',SNRs(:,3),'g')
subplot(2,1,2);plot(SNRms(:,1),'r',SNRms(:,2),'b',SNRms(:,3),'g')
end

