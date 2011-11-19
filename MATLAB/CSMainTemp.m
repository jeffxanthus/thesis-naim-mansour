function [result] = CSMainTemp(signal, fs)
%CSMain - Main function of the CS Declipping program.
%   signal: The input audio signal.
%   fs:     Sample rate in Hz. If not specified, 44.1kHz is chosen.
%Naim Mansour

%If no sample frequency specified
tic;
if nargin<2
    fs=44.1*10^3;
end
%Some checks to make sure the input has the correct dimensions.
[rs cs]=size(signal);
if ~(rs==1 || cs==1)
    disp('Input signal needs to be a vector.')
    return;
end
if rs<cs
    signal=signal';
end
[rs cs]=size(signal);

%Parameter selection
frameLength=20; %in milliseconds
if mod(fs*frameLength,1000)~=0
    disp('Infeasible sampling frequency, using default')
    fs=44.1*10^3;
end
frameOverlap=50; %percentage -- FOR NOW, USE ONLY 50%, OTHERWISE BARTLETT DOESN'T WORK
numberOfSamples=fs*frameLength/1000; %also minimal signal length

%Rough approach, to be changed later
if rs<numberOfSamples
    disp('Signal too short')
    return;
end

%To avoid rounding issues, and synthesis frame issues
if mod(numberOfSamples*(frameOverlap/100),1)~=0
    disp('Illegal overlap ratio')
    return;
end

%Keep only the part of the signal divisible into subframes
signal(end-(mod(rs,numberOfSamples)-1):end,:)=[];
[rs cs]=size(signal);

%Division into frames
disp('Dividing into frames...')
T=[];
position=1;
i=1;
while position+(numberOfSamples)<=rs
    T(:,i)=signal(position:position+(numberOfSamples),1);
    position=position+(numberOfSamples)-((frameOverlap/100)*numberOfSamples);
    i=i+1;
end


U=[];
[rst cst]=size(T);
% Could be combined with previous for-loop
% Apply the declipping to each frame
disp('Declipping...')
% for j=1:cst
%     disp(['Now declipping frame ' int2str(j) ' out of ' int2str(cst)])
%     U(:,j)=CSDeclip(T(:,j));
%     subplot(2,1,1);plot(U(:,j),'-');
%     subplot(2,1,2);plot(T(:,j),'-');
% end
U=T;
disp('Reconstructing...')
%Reconstruction (should work for any overlap ratio and compatible window)
window=hann(numberOfSamples+1); %Hann window (with 50% overlap)
shiftAmount=(1-(frameOverlap/100))*numberOfSamples; %=440 in typical case
%Condition on first frame
result=[U(1:shiftAmount,1)' (window(shiftAmount+1:end,1).*U(shiftAmount+1:end,1))'];
position=1;
for k=2:cst
    if k==cst
        %Condition on last frame
        currentBlock=[(window(1:shiftAmount,1).*U(1:shiftAmount,k))' U(shiftAmount+1:end,k)'];
    else
        currentBlock=(window.*U(:,k))';
    end
    unalteredPart=result(1,1:position+(shiftAmount-1));
    summedPart=result(1,position+shiftAmount:end)+currentBlock(1,1:end-shiftAmount);
    unalteredBlockPart=currentBlock(1,end-(shiftAmount-1):end);
    
    result=[unalteredPart summedPart unalteredBlockPart];
    position=(k-1)*shiftAmount+1;
end
toc;

size(result)
size(signal)
subplot(2,1,1);plot(signal);
title('Clipped signal')
subplot(2,1,2);plot(result-signal);
title('Reconstructed signal')
%Mystical error when no declipping - no perfect frame reconstruction...
end
