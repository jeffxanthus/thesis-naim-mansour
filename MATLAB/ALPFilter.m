function [newReconstruction] = ALPFilter(input, reconstruction)
%ALPFILTER Summary of this function goes here
%   Detailed explanation goes here

[a b]=size(input);
[c d]=size(reconstruction);
if a~=1
    input=input';
end
if c~=1
    reconstruction=reconstruction';
end

origSamples=[];
MaxS=max(input);
MinS=min(input);
ratio=determineLPFactor(reconstruction);
ratio

for t=1:length(input) 
    if(input(1,t)>=MaxS || input(1,t)<=MinS)...
        && (((t==1 && input(1,t+1)-input(1,t)==0)  || (t==length(input) && input(1,t)-input(1,t-1)==0))...  
        || ((t~=1 && t~=length(input)) && input(1,t-1)-input(1,t+1))==0)
        origSamples=[origSamples t];
    end
end

x=dct(reconstruction, length(reconstruction));
cutOff=round(ratio*length(reconstruction));
xNew=x(1,1:cutOff);
reconstruction=idct(xNew,length(reconstruction));

% LP only on clipped samples
% reconstruction(:,origSamples)=reconstruction(:,origSamples);
newReconstruction=reconstruction;
end

