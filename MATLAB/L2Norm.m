function[f]=L2Norm(x) 
%L2NORM Takes the L2 norm of the given vector
global B
global yU

tol=0.1;
[a b]=size(x);
if((a~=1)&&(b~=1))
    disp('Input is not a vector')
    return;
end
if(a==1)
    x=x';
end

f=norm(yU-B*x);
end

