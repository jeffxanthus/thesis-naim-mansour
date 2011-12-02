function [output] = L0Norm(data)
%L0NORM Takes the L0-norm of a given vector

[x y]=size(data);
if((x~=1)&&(y~=1))
    disp('Input is not a vector')
    return;
end
if(x==1)
    data=data';
end

output=0;

for i=1:length(data)
    if(data(i,1)>0)
        output=output+1;
    end
end
end