function [dataOrig1, reconstruction1, reconstruction2] = AxBeComparisonTest(clip,size)
%AXBECOMPARISONTEST 
% Author: Naim Mansour
global methodChoice

methodChoice=3;
offSet=(round(rand()*20000))+(round(rand()*100000));

[dataOrig1, reconstruction1, origSamples1, missingSamples1, SNR1, SNRm1, SNROrig]=Simulation(clip,size,1,offSet);
[dataOrig2, reconstruction2, origSamples2, missingSamples2, SNR2, SNRm2, SNROrig]=Simulation(clip,size,2,offSet);

subplot(4,2,1);plot(origSamples1,'r.')
title('Clipped values in original signal (Ax=y)')
subplot(4,2,2);plot(origSamples2,'r.')
title('Clipped values in original signal (Ax+Be=y)')
subplot(4,2,3);plot(missingSamples1,'.')
title('"Clipped values" in reconstructed signal (Ax=y)')
subplot(4,2,4);plot(missingSamples2,'.')
title('"Clipped values" in reconstructed signal (Ax+Be=y)')
subplot(4,2,5);plot(abs((origSamples1-missingSamples1)./origSamples1))
title('Relative error for Ax=y')
axis([0 length(origSamples2) 0 0.5])
subplot(4,2,6);plot(abs((origSamples2-missingSamples2)./origSamples2))
title('Relative error for Ax+Be=y')
axis([0 length(origSamples2) 0 0.5])
subplot(4,2,7);plot(reconstruction1)
title('Reconstruction with Ax=y')

subplot(4,2,8);plot(reconstruction2)
title('Reconstruction with Ax+Be=y')

SNROrig

SNR1
SNR2

SNRm1
SNRm2
end

