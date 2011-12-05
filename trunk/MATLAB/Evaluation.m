function [SNR,ODG] = Evaluation(signal, reconstruction, fs, noBits)
%EVALUATION 
% SNR evaluation of the reconstruction.
% Author: Naim Mansour
addpath('C:\Users\Naim\Documents\K.U. Leuven\Thesis\MATLAB\PEAQ');

if nargin<4
    noBits=16;
end
SNR=10*log10(norm(signal,2).^2/norm(signal-reconstruction,2).^2);

wavwrite(reconstruction,48000,noBits,'reconstruction.wav');
wavwrite(signal,48000,noBits,'signal.wav');
ODG=PQEvalAudio('signal.wav','reconstruction.wav')
end

