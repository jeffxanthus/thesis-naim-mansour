function [SNR] = perceptualTest(signal, method, clippinglevel)
% testding voor de perceptuele dingen
% Steven De Hertogh

global methodChoice 

methodChoice = 3;

sig = wavread(signal);
sig = sig(:,1);
sig = sig(20:6000);
result = CSPerceptualMain(Clipold(sig,clippinglevel),method, 44000);
SNR = Evaluation(sig', result, 44100, 16);

end