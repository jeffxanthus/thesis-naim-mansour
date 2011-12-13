function A = spdiags(B,d,m,n)
% tomSym/spdiags - Overloaded function
%
% A = SPDIAGS(B,d) extracts diagonals of B specified by d.
% A = SPDIAGS(B,d,A) modifies a sparse matrix.
% A = SPDIAGS(B,d,m,n) creates an m-by-n sparse matrix.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargout>1
    error('This calling syntax for spdiags is not supported by tomSym.');
end

if nargin==2
    % A = SPDIAGS(B,d) - Extract diagonals of B specified by d.
    Ac = cell(1,length(d));
    n = min(size(B));
    for i=1:length(d)
        if d(i)==0 && size(B,1)==size(B,2)
            Ac{i} = getdiag(B);
        else
            Ac{i} = vec(lookup(B,find(spdiags(ones(n,1),d(i),size(B,1),size(B,2)))));
            if d(i)>0 || size(B,1)~=size(B,2) % Stupid Matlab bug
                Ac{i} = [zeros(n-length(Ac{i}),1); Ac{i}];
            elseif d(i)<0
                Ac{i} = [Ac{i}; zeros(n-length(Ac{i}),1)];
            end
        end
    end
    A = horzcat(Ac{:});
    return
end

if nargin>=3
    % A = SPDIAGS(B,d,A) modifies a sparse matrix.
    % A = SPDIAGS(B,d,m,n) creates an m-by-n sparse matrix.
    
    if nargin==4
        A = spzeros(m,n);
    end
    
    n = min(size(A));
    
    for i=1:length(d)
        if d(i)==0 && size(A,1)==size(A,2)
            A = A + setdiag(submatrix(B,1:size(A,1),i));
        else
            if d(i)>0 || size(B,1)~=size(B,2) % Stupid Matlab bug
                [ii,jj] = find(spdiags(ones(n,1),d(i),size(A,1),size(A,2)));
                A = A + sparse(ii,jj,submatrix(B,n-length(ii)+1:n,i),size(A,1),size(A,2));
            elseif d(i)<0
                [ii,jj] = find(spdiags(ones(n,1),d(i),size(A,1),size(A,2)));
                A = A + sparse(ii,jj,submatrix(B,1:length(ii),i),size(A,1),size(A,2));
            end
        end
    end
    return
end
   
        
