function [result,missingSamples] = CSPerceptualMain(signal, method, fs)
%CSMain - Main function of the CS Declipping program.
%   signal: The input audio signal.
%   fs:     Sample rate in Hz. If not specified, 44.1kHz is chosen.
%   method: Type 1 for classical y=Ax constraints, 2 for y=Ax+Beconstraints
%                3 for perceptual y =Ax
%   Naim Mansour and Steven De Hertogh

%If no sample frequency specified
tic;
if nargin<3
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
fs=44000; %Artificial - testing purposes

%%NEW APPROACH - ADAPTIVE FRAMELENGTH - CODE BULKY SO FAR, to be altered
%%later
clipped1=length(find(signal == max(signal)));
clipped2=length(find(signal == min(signal)));

maxAmountOfFrames=0;
if ~(clipped1==1 && clipped2==1)
    maxAmountOfFrames=floor((rs-(clipped1+clipped2))/(2*estimateSparsity(signal)-1).^2)
end


%%TODO: choose frame length according to constraint, but also to
%%numerical (UPPER)/psychoaccoustic (LOWER) constraints  --explain in text

%Parameter selection - to be discretised later, fs only necessary for
%playback and PEAQ
frameLength=floor(((rs/maxAmountOfFrames)/fs)*1000)
%Computational bound
frameLength=min(frameLength,60)

% frameLength=30; %in milliseconds
if mod(fs*frameLength,1000)~=0
    disp('Infeasible sampling frequency, using default')
    fs=44.1*10^3;
end
frameOverlap=50; %percentage -- FOR NOW, USE ONLY 50%, OTHERWISE BARTLETT DOESN'T WORK
numberOfSamples=fs*frameLength/1000 %also minimal signal length

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
origSig = signal;
signal(end-(mod(rs,numberOfSamples)-2):end,:)=[];
[rs cs]=size(signal);


%Division into frames
disp('Dividing into frames...')
T=[];
U=[];

position=(1-(frameOverlap/100))*numberOfSamples+1;
window=hann(numberOfSamples+1); %Hann window (with 50% overlap)
shiftAmount=(1-(frameOverlap/100))*numberOfSamples+1; %=442 in typical case
masking = maskingThreshold(1);
i=1;
while position+(numberOfSamples)<=rs
    if i == 1
        T(:,i)=signal(1:numberOfSamples+1,1);
    else
        T(:,i)=signal(position:position+numberOfSamples,1);
    end
    disp(['Now declipping frame ' int2str(i) ])
    if(method==1)
        [dummy U(:,i)]=CSDeclip(T(:,i));
    elseif (method==2)
        U(:,i)=CSDeclipAlternate(T(:,i));
    elseif (method==3)
        [dummy U(:,i)]=CSPerceptualDeclip(T(:,i), masking);
    end
    
    disp('Reconstructing...')    
    %Condition on first frame
    if i == 1
        %if cst==1 %If only 1 frame
        %    result=U';
        %else
            result=[U(1:shiftAmount-1,1)' (window(shiftAmount:end,1).*U(shiftAmount:end,1))'];
       % end
    end
    if i==rs
        %Condition on last frame
        currentBlock=[(window(1:shiftAmount,1).*U(1:shiftAmount,i))' U(shiftAmount+1:end,i)'];
    else
        currentBlock=(window.*U(:,i))';
    end
    unalteredPart=result(1,1:position-1);
    summedPart=result(1,position:end)+currentBlock(1,1:shiftAmount);
    unalteredBlockPart=currentBlock(1,shiftAmount+1:end);
    result=[unalteredPart summedPart unalteredBlockPart];
    masking = meanMaskingThreshold(result);
    position=position+shiftAmount-1;
    i=i+1;
end



%Non-multiple part
if(method==1)
        [dummy nonMultipleRec]=CSDeclip(nonMultiplePart);
    else
        nonMultipleRec=CSDeclipAlternate(nonMultiplePart);
end

result=[result nonMultipleRec'];
toc;

missingSamples=[];
for u=Miss
    try
        missingSamples=[missingSamples result(1,u)];
    catch ME
        disp('Missing sample skipped');
    end
end

maskingThreshold(result);

figure();
subplot(3,1,1);plot(origSig);
title('Clipped signal')
subplot(3,1,2);plot(result);
title('Reconstructed signal')
% subplot(3,1,3);plot(abs(signal-result));
subplot(3,1,3);plot(missingSamples);
title('Reconstructed missing samples');
end