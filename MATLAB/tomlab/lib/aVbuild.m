% function alphaV=aVbuild(A,Its)
%
% Build matrix alphaV with line search steps, from
% matrix A with rows: iterations step # line search length
% Its is total number of iterations
% Default Its is the maximal iteration number in matrix A

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1997-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Nov 6, 1997.   Last modified Jan 10, 2002.

function alphaV=aVbuild(A,Its)

if isempty(A) 
   alphaV=[];
   return
end
mx=max(A(:,1));

if isempty(mx), mx=0; end

if nargin < 2, Its=[]; end
if isempty(Its), Its=mx; end

m=size(A,1);

if m==0 & size(A,2)==2
   % All steps were unit steps
   alphaV=ones(Its,1);
elseif mx==0
   alphaV=ones(Its,1);
   fprintf('ERROR in aVbuild: size(A)=%d %d; Its %d\n',size(A),Its);
else
   maxa=max(sum([ ones(m,1)*[1:mx]==A(:,1)*ones(1,mx) ]));
   alphaV=[ones(Its,1),zeros(Its,maxa)];
   old=0; j=0;
   for i=1:m
       if A(i,1)==old
          j=j+1;
       else
          j=2;
          old=round(real(A(i,1)));
       end
       alphaV(old,j)=A(i,2);
   end
end