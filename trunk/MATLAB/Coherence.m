function [coherence, coherence_un, coherence_sq] = Coherence(A,B)
%COHERENCE
% coherence: unnormalized coherence
% coherence_n: normalized coherence through division by product of 2-norms
% coherence_sq: unnormalized coherence multiplied with sqrt(n)
%Naim Mansour
coherence=0;
coherence_un=0;
coherence_sq=0;
index=[1 1];
%Coherence of matrix A
if nargin < 2,
    disp('a')
    [m n]=size(A);
    for i=1:n-1
        for j=i+1:n
            newcoh=(abs(dot(A(:,i),A(:,j))));
            if newcoh>coherence
                coherence=newcoh;
                if coherence<1e-10 coherence=0;end
                index=[i j];
            end
        end
    end
    k=index(1,1);
    l=index(1,2);
    coherence_n=coherence/(norm(A(:,k),2)*norm(A(:,l),2));
    coherence_sq=coherence*sqrt(n);
    
elseif nargin >= 2,
    [m1 n1]=size(A);
    [m2 n2]=size(B);
    for i=1:n1
        for j=1:n2
            newcoh=(abs(dot(A(:,i),B(:,j))));
            if newcoh>coherence
                coherence=newcoh;
                if coherence<1e-10 coherence=0;end
                index=[i j];
            end
        end
    end
    k=index(1,1);
    l=index(1,2);
    coherence_n=coherence/(norm(A(:,k),2)*norm(B(:,l),2));
    coherence_sq=coherence*sqrt(n1);
end

disp('Maximal recoverable sparsity level, less than')
0.5*(1+1/coherence_n)
disp('If unitary, satisfies RIP of order k with delta=(k-1)*mu for all')
disp('k<');1/coherence_n

end

