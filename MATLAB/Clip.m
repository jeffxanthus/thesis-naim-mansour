function [output] = Clip(input,threshold,clippedAmount)
%CLIP - Naim Mansour
%threshold el. of [0,1]
%Naim Mansour
limit=threshold*max(input);
[rs cs]=size(input);
x=input;
    
if nargin == 2 
    if(rs>1)
        input=input';
    end
    for i=1:length(input)
        if(input(1,i)>limit)
            input(1,i)=limit;
        end
         if(input(1,i)<-limit)
            input(1,i)=-limit;
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

subplot(2,1,1);plot(x,'.');
axis([0 length(x) min(x)-1 max(x)+1])
subplot(2,1,2);plot(output,'.');
axis([0 length(x) min(x)-1 max(x)+1])
end

