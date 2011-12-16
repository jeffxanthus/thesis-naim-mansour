function [data,largeData,mediumData,smallData,tinyData,fs,noBits] = InitializeTestVariables(fileName,offSet)
%INITIALIZETESTVARIABLES 
% Author: Naim Mansour
addpath('../Samples')

[data,fs,noBits]=wavread(fileName);

[rs cs]=size(data);
if rs~=1
    if rs==2
        data=data(:,1)';
    else
        data=data';
    end
end

% i=1;
% remove=[];
% while data(1,i)==0
%     remove=[remove i];
%     i=i+1;
% end
% remove
% data(1,remove)=[];

data=data(1,20000:end);

%Temporary solution, only testing purposes
largeData=[];%data(1,offSet+10001:offSet+2000000); %1990000 samples
mediumData=data(1,offSet+10001:offSet+510000); %500000 samples
smallData=data(1,offSet+20001:offSet+30000); %10000 samples
tinyData=data(1,offSet+10001:offSet+13000);  %3000 samples
end

