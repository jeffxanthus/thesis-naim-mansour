function [SNR SNRCLIP] = optimizationTest(signal)
% testding voor de perceptuele dingen
% Steven De Hertogh
global methodChoice

methodChoice = 6;

sig = wavread(signal);
sig = sig(:,1);
sig = sig(200:800);
result = CSPerceptualDeclip(Clipold(sig, 0.7)');
SNR=10*log10(norm(sig,2).^2/norm(sig-idct(result),2).^2)
SNRCLIP = 10*log10(norm(sig,2).^2/norm(sig-Clipold(sig,0.7)',2).^2)


end