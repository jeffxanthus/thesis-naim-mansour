function [clips] = ClipTest(fileName)
%CLIPTEST 

for t=1:5
    [data,largeData,mediumData,smallData,tinyData,fs,noBits] = InitializeTestVariables(fileName,t*100000);

    data=smallData;
    clipVec=[0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2];
    originalSize=length(data);
    newSize=zeros(8,1);

    for i=1:8
        clippedData=Clip(data,clipVec(1,i));
        for j=1:originalSize
            if data(1,j)-clippedData(1,j)~=0
                newSize(i,1)=newSize(i,1)+1;
            end
        end
        newSize(i,1)=newSize(i,1)./originalSize;
    end

    clips=[clipVec ; newSize']
    pause
end
