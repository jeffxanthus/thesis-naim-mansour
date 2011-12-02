function [output] = L1Norm(data)
%L1NORM Takes the L1-norm of a given vector

[x y]=size(data);
if((x~=1)&&(y~=1))
    disp('Input is not a vector')
    return;
end
if(x==1)
    data=data';
end
output=norm(data,1);
end

