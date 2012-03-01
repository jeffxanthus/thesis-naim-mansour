function [SNR SNRCLIP] = optimizationTest(signal, optim)
% testding voor de perceptuele dingen
% Steven De Hertogh
% global methodChoice
% 
% methodChoice = 6;

sig = wavread(signal);
sig = sig(:,1);
sig = sig(200:800);
clipped = Clipold(sig, 0.5)';
tic
result = CSPerceptualDeclip(clipped, maskingThreshold(sig'), optim);
toc

size(sig)
size(result)

SNR=10*log10(norm(sig,2).^2/norm(sig-idct(result),2).^2)
SNRCLIP = 10*log10(norm(sig,2).^2/norm(sig-clipped,2).^2)

result = idct(result);

figure();
subplot(3,1,1);plot(clipped);
title('Clipped signal')
subplot(3,1,2);plot(result);
title('Reconstructed signal')
subplot(3,1,3);plot(sig);
title('Original signal')

end