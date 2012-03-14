function [output] = Clip(input,threshold)
%CLIP - Naim Mansour
%threshold el. of [0,1]
%Naim Mansour

%Avoid outliers messing up the clipping (misleading results)
%MAX case
newMax=max(input);
newMin=min(input);

% [sortedSig indicesSig]=sort(input,'descend');
% meanSig=mean(sortedSig);
% stdSig=std(sortedSig)
% t=1;
% %Could be more detailed: later
% while (sortedSig(1,t)>meanSig+3*stdSig)
%     t=t+1;
% end
% newMax=sortedSig(1,t);
% 
% t=length(sortedSig);
% while (abs(sortedSig(1,t))>meanSig+3*stdSig)
%     t=t-1;
% end
% newMin=sortedSig(1,t);

% limit=threshold*max(input);
upperLimit=threshold*newMax;
lowerLimit=threshold*newMin;
[rs cs]=size(input);
x=input;

if nargin == 2 
    if(rs>1)
        input=input';
    end
    for i=1:length(input)
        if(input(1,i)>upperLimit)
            input(1,i)=upperLimit;
        end
         if(input(1,i)<lowerLimit)
            input(1,i)=lowerLimit;
        end
    end
end

if nargin > 2
    if mod(clippedAmount,2)~=0
        clippedAmount=clippedAmount+1;
    end
    [rs cs]=size(input);
    if rs==1
        input=input';
    end
    [inputSorted indices]=sort(input,'descend'); %dubious, only for test purposes
    size(input)
    size(inputSorted)
    %Upper bound
    toReplace=indices(1:clippedAmount/2,1)';
    for t=toReplace
        input(t,1)=inputSorted(clippedAmount/2,1); 
    end
    %Lower bound
    toReplace2=indices(end-(clippedAmount/2)+1:end,1)';
    for v=toReplace2
        input(v,1)=inputSorted(end-(clippedAmount/2)+1,1); 
    end
end
output=input;

% subplot(2,1,1);plot(x,'.');
% axis([0 length(x) min(x)-1 max(x)+1])
% subplot(2,1,2);plot(output,'.');
% axis([0 length(x) min(x)-1 max(x)+1])
end

