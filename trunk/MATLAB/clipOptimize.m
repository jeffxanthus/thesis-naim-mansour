function [clip, signal] = clipOptimize(fileName,PEAQDesired)
%CLIPOPTIMIZE Summary of this function goes here
%   Detailed explanation goes here
fs=44100;
inputSignal=wavread(fileName);
PEAQ=5;
i=0.99858;
rateOfChange=0.0000001;
while PEAQ>PEAQDesired
    if(PEAQ-PEAQDesired<rateOfChange)
        rateOfChange=0.01*rateOfChange;
    end
    i=i-rateOfChange;
    clippedSignal=Clip(inputSignal,i)';
%     sound(clippedSignal,44100);
%     spause
    [dummy PEAQ]=Evaluation(inputSignal, clippedSignal,fs,16);
    PEAQ
    i
end

end

