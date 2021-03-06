function [result,missingSamples] = CSMain(signal, method, fs)
%CSMain - Main function of the CS Declipping program.
%   signal: The input audio signal.
%   fs:     Sample rate in Hz. 
%   method: Type 1 for classical y=Ax constraints, 2 for y=Ax+Beconstraints
%Naim Mansour
% addpath('../')
%If no sample frequency specified
global methodChoice
global clip
global fL
global LPFactor
tic;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OLD CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fs=44000; %Artificial - testing purposes
%%NEW APPROACH - ADAPTIVE FRAMELENGTH - CODE BULKY SO FAR, to be altered
%%later
% clipped1=length(find(signal == max(signal)));
% clipped2=length(find(signal == min(signal)));
% maxAmountOfFrames=0;
% if ~(clipped1==1 && clipped2==1)
%     maxAmountOfFrames=floor((rs-(clipped1+clipped2))/(2*estimateSparsity(signal)-1).^2);
% end
%Parameter selection - to be discretised later, fs only necessary for
%playback and PEAQ
% frameLength=floor(((rs/maxAmountOfFrames)/fs)*1000);
%Computational bound
% frameLength=min(frameLength,60); %Problems with divisibility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Make sure fs equals 44kHz.
if fs~=44000
    disp('Sampling frequency has to be 44kHz.')
    return;
end

%No frame lengths under 30ms
frameLength=max(fL,30);

if mod(fs*frameLength,1000)~=0
    disp('Infeasible sampling frequency')
    return;
end
frameOverlap=50; %percentage --
numberOfSamples=fs*frameLength/1000; %also minimal signal length - 1323

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

%Used to be able to calculate SNRm
Miss=[];
MaxS=max(signal);
MinS=min(signal);
for t=1:rs %border values possibly added to the clipped values -OK NOW
    if (signal(t,1)>=MaxS || signal(t,1)<=MinS)...
        && (((t==1 && signal(t+1,1)-signal(t,1)==0)  || (t==rs && signal(t,1)-signal(t-1,1)==0))... 
        || ((t~=1 && t~=rs) && signal(t-1,1)-signal(t+1,1))==0)
        Miss=[Miss t];
    end
end

%Keep only the part of the signal divisible into subframes, the other part
%will be processed separately (temp, maybe better approach later)
nonMultiplePart=signal(end-(mod(rs,numberOfSamples)-2):end,:);
origSignal=signal;
signal(end-(mod(rs,numberOfSamples)-2):end,:)=[];
[rs cs]=size(signal);


%Division into frames
disp('Dividing into frames...')
T=[];
T(:,1)=signal(1:numberOfSamples+1,1);
position=(1-(frameOverlap/100))*numberOfSamples+1;
i=2;
while position+(numberOfSamples)<=rs
    T(:,i)=signal(position:position+numberOfSamples,1);
    position=position+(1-(frameOverlap/100))*numberOfSamples;
    i=i+1;
end


U=[];
[rst cst]=size(T);
% Could be combined with previous for-loop

% Apply the declipping to each frame
disp('Declipping...')
for j=1:cst
    disp(['Now declipping frame ' int2str(j) ' out of ' int2str(cst)])
    if(method==1)
        [dummy U(:,j)]=CSDeclip(T(:,j));
    else
        [U(:,j) dummy]=CSDeclipAlternate(T(:,j));
    end
    subplot(2,1,1);plot(U(:,j),'-');
    subplot(2,1,2);plot(T(:,j),'-');
end
%Non-multiple part
if(method==1)
        [dummy nonMultipleRec]=CSDeclip(nonMultiplePart);
    else
        [nonMultipleRec dummy]=CSDeclipAlternate(nonMultiplePart);
end
%   U=T;
%   nonMultipleRec=nonMultiplePart; 
disp('Reconstructing...')
%Reconstruction (should work for any overlap ratio and compatible window)
window=hann(numberOfSamples+1); %Hann window (with 50% overlap)
shiftAmount=(1-(frameOverlap/100))*numberOfSamples+1; %=442 in typical case
%If only 1 frame
if cst==1
    result=U';
else
    %Condition on first frame
    result=[U(1:shiftAmount-1,1)' (window(shiftAmount:end,1).*U(shiftAmount:end,1))'];
end
position=shiftAmount;
for k=2:cst
    if k==cst
        %Condition on last frame
        currentBlock=[(window(1:shiftAmount,1).*U(1:shiftAmount,k))' U(shiftAmount+1:end,k)'];
    else
        currentBlock=(window.*U(:,k))';
    end
    unalteredPart=result(1,1:position-1);
    summedPart=result(1,position:end)+currentBlock(1,1:shiftAmount);
    unalteredBlockPart=currentBlock(1,shiftAmount+1:end);
    result=[unalteredPart summedPart unalteredBlockPart];
    position=position+shiftAmount-1;
end

result=[result nonMultipleRec];

if LPFactor
    result=ALPFilter(signal,result);
end
toc;

missingSamples=[];
for u=Miss
    try
        missingSamples=[missingSamples result(1,u)];
    catch ME
        disp('Missing sample skipped');
    end
end

subplot(3,1,1);plot(origSignal);
title('Clipped signal')
subplot(3,1,2);plot(result);
title('Reconstructed signal')
subplot(3,1,3);plot(missingSamples);
title('Reconstructed missing samples');
end