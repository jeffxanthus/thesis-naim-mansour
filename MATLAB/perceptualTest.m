function [SNR] = perceptualTest(signal, method, clippinglevel, optim)
% testding voor de perceptuele dingen
% Steven De Hertogh

global methodChoice 

methodChoice = method;

sig = wavread(signal);
sig = sig(:,1);
sig = sig(20:6000);
result = CSPerceptualMain(Clip(sig,clippinglevel),method, 44000, optim);
SNR = Evaluation(sig', result, 44100, 16);

end