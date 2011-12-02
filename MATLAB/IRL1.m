function [x] = IRL1(A, y, N, maxIters, lambda, OptTol)
%IRL1  - Iteratively reweighted L1  (development in progress)
%Naim Mansour
global fFf

iterCount=3; %If this is set to 1, unweighted L1 minimization is carried out.
[rs cs]=size(A);
%Initialization
W=diag(ones(N,1));
it=0;
xNew=1*10^5.*ones(N,1);
eps=1*10^(-4);

x=zeros(N,1);
% x0=SolveBP(A,y,N,maxIters,lambda,OptTol);
options = struct('Verbose', 0);
[x0,niter,residuals,outputData,opts]=NESTA(A,[],y,0.02,1e-3,options);
fFf=0.95;
epsilon=chooseEpsilons(x0,N,fFf);

%If BP: min lambda*||x||1 + ||Ax-b||2
%Weighted BP: y=Wx, x=W^(-1)x, min lambda*||y||1 + ||AW^(-1)y-b||2
while (it<iterCount) %norm(x-xNew)>eps &&
    disp(['Iteration count is now ' int2str(it)])
    x=xNew;
%     B=A/W;
%     xNew=SolveBP(B,y,N,maxIters,lambda,OptTol);
    %GET BACK TO 5dB INCREASE PERFORMANCE!!!
    options = struct('Verbose', 0, 'U', W,'maxiter',250,'MaxIntIter',3,'xplug',x0);
    [xNew,niter,residuals,outputData,opts]=NESTA(A,[],y,0.02,1e-3,options);
%     xNew=W\xNew;
    close all;
%     subplot(2,1,1);plot(xNew)
%     subplot(2,1,2);plot(idct(xNew))
%     pause
    xNew=xNew';
    newWeights=(abs(xNew)+epsilon').^(-1);
    W=diag(newWeights);
    it=it+1;
end
x=xNew';

    %Not used yet
    function[epsilon]=chooseNewEpsilon(x)
        %- Chooses new epsilon based on method proposed by Candès et al.
        xSorted=sort(abs(x),2,'descend');

        i0=round(rs/(4*log(N/rs))); %round introduces because integer indexing - effect?
        xBound=x(1,i0); %Index too high, since rs is almost equal to N
        epsilon=max(xBound,10^(-3));
    end

    function[epsilon]=chooseEpsilons(x0,N,f)
        epsilon=[];
        fFf=0.95;
        for i=1:N
              epsilon(i,1)=(f)*x0(i,1);
        end
    end
end

