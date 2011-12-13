function H = prodJ1J1(M,dim)
% tomSym/prodJ1J1 - second derviative of prod
%
% J = prodJ1J1(M,dim) computes the second derivative matrix of the
% function call prod(M,dim)
%
% Each row of the returned matrix J corresponds to an element of
% prodJ1(M,dim) and each column corresponds to an element of M.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if dim==2
    M = M';
end

p = prod(M);

[m,n] = size(M);

H = spalloc(n*(n*m),n*m,(n*m)*(m-1));

for i=1:n
    if p(i)==0
        % At least one zero among the factors.
        j = find(M(:,i)==0);
        if length(j)>2
            % All derivs = 0
            continue;
        elseif length(j)==2
            v = p(:,i);
            v(j) = 1;
            v = prod(v);
            HM = zeros(m,m);
            HM(j(1),j(2)) = v;
            HM(j(2),j(1)) = v;
        else %length(j)==1
            v = M(:,i);
            v(j) = 1;
            v = prod(v)./v;
            HM = zeros(m,m);
            HM(:,j) = v;
            HM(j,:) = v';
            HM(j,j) = 0;
        end
    else
        % No zeros among the factors
        HM = repmat(p(i)./M(:,i),1,m)./repmat(M(:,i)',m,1);
        HM(sub2ind([m,m],1:m,1:m)) = 0;
    end

    if dim==1
        ix = i+(i-1)*m*n+(0:(m-1))*n;
        jx = (i-1)*m+(1:m);
    else
        ix = i+(i-1)*n+(0:m-1)*n*n;
        jx = i+(0:m-1)*n;
    end
    
    H(ix,jx) = HM;
end

