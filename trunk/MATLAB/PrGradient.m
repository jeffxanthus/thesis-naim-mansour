function [xSol] = PrGradient(A,y,B,theta,x0,maxIter)
%PRGRADIENT 
%M=number of samples, N=total signal length, L=missing sample length
%A: DCT base (size MxN)
%y: samples (size Mx1)
%B: Constraint matrix (size LxN)
%theta: clipping constraints (size Lx1)
%x0: Initial x (size Nx1)
%maxIter: maximum amount of iterations

i=1;
x=1e5;
xNew=x0;
eps=0.01;   %accuracy
alpha=0.01; %step-size

while i<maxIter && norm(xNew-x)>eps
    x=xNew;
    C=(B*x-theta);
    if C<=zeros(length(theta),1)
        g=gradient(norm(y-A*x).^2);
        %g will have size Mx1 in this case
    else
        g=gradient(C);
        %g will have size Lx1 in this case
    end     
    xNew=x-alpha.*g;
    %g needs to have size Nx1 in order to be able to do this... ??
end

xSol=x;
end

