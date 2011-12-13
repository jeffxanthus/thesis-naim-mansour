% ====================================================================
function X = randEllips(x_L,x_U,proc,n,dist,X)
% ====================================================================

d = length(x_L);
if nargin < 5
   dist = [];
   if nargin < 4
      n = 2*d;
   end
end
xD = x_U-x_L;
if isempty(dist), dist=norm(xD); end

minDist=proc/100*dist;

%X=rand(d,1).*dist+x_L;

if nargin < 6
   X            = [x_L + xD/2,zeros(d,n-1)];   % center point
   k            = 1; 
else
   k            = size(X,2);
   X(:,k+1:k+n) = 0;   % Initialize
   n            = n+1; % Increase n, because we are looping to n-1
end

xNew = x_L + rand(d,1).*xD;

for i=1:n-1
    j = 0;
    %[t,m]=size(X);
    %while min(sqrt(sum((xNew*ones(1,m)-X).^2)))<minDist & j < 1000
    while min(tomsol(30,xNew,X)) < minDist & j < 1000
        xNew = x_L + rand(d,1).*xD;
        j    = j+1;
    end
    X(:,k+i) = xNew;
end