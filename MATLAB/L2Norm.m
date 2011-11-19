function[c,ceq]=L2Norm(x) 
%L2NORM Takes the L2 norm of the given vector
global A
global samples

tol=0.01;
[a b]=size(x);
if((a~=1)&&(b~=1))
    disp('Input is not a vector')
    return;
end
if(a==1)
    x=x';
end

c=norm(A*x-samples)-tol;
ceq=[]; 

end

