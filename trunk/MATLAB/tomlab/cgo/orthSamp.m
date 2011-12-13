% function X = orthSamp(n, s, x_L, x_U, method, PriLev)
%
% Orthogonal Sampling:  In each dimension, the sample space is divided into s 
%                       equally probable subspaces. Every subspace block s^d
%                       will have n points in it. In total, N = n*s^d points.
%
%
%       N=n*s^d   OS:   |---|---|---|---|          LH:   |---|---|---|---|
%                       |   | X |   |   |                | X |   |   |   |
%       N=4             |---|---|---|---|                |---|---|---|---|
%                       |   |   | X |   |                |   | X |   |   |
%       n=1             |---|---|---|---|                |---|---|---|---|
%       s=2             | X |   |   |   |                |   |   |   | X |
%       d=2             |---|---|---|---|                |---|---|---|---|
%                       |   |   |   | X |                |   |   | X |   |
%                       |---|---|---|---|                |---|---|---|---|


% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written April 10, 2008.    Last modified Aug 24, 2008.

function X = orthSamp(n, s, x_L, x_U, method, PriLev)

d  = length(x_L);

if method < 0
   IntVals = 1;
   method = abs(method);
else
   IntVals = 0;
   xD = x_U-x_L;
end


if isempty(s)
   N = n;      % Interpret input as the desired number of sample points.
   % Find suitable values for n and s.
   n = 1:N-1;
   s = ceil((N./n).^(1/d));
   opt = n.*s.^d;
   [v,I] = min(opt-N);
   n = n(I);
   s = s(I);
   if v ~= 0
      if PriLev > 2
         fprintf(' Not possible to sample %d points in %d dimensions.\n',N,d)
      end
      N = n*s^d;
      if PriLev > 2
         fprintf(' New number of sample points %d.\n',N)
      end
   end
end

% Useful parameters
D = s^(d-1);   % D = s^(d-1)
S = s*D;       % S = s^d
M = n*D;       % M = n*s^(d-1);
N = n*S;       % N = n*s^d;

X = zeros(d,N);

if PriLev > 1
   fprintf(' Orthogonal Sampling: %d points  ( n=%d, s=%d, d=%d ).\n',N,n,s,d)
end

if method == 1;

   A = repmat(1:N,d,1);
   x = zeros(1,d);
   
   BOX = zeros(N,d);
   
   p = N;
   q = 1;

   for k = 1:d
      p = p/s;
      q = q*s;
      kk = 1:q;
      tmp = repmat(mod(kk-1,s)+1,p,1);
      BOX(:,k) = tmp(:);
   end
   
   K = randperm(N);
   for k = K
      for kk = 1:d
         boxNr = BOX(k,kk);
         alt = intersect( (boxNr-1)*M+1:boxNr*M , A(kk,:) );
         % pick first alternative
%          x(kk) = alt(1);
         % randomly pick one of the alternatives
         L = length(alt);
         I = round((L-1)*rand(1)+1);
         x(kk) = alt(I);
         % remove alternative from A
         %I = find(A(kk,:) == x(kk));
         %A(kk,I) = 0;
         A(kk, A(kk,:) == x(kk)) = 0;
      end
      X(:,k) = x;
      
      A = sort(A,2);
      A(:,1) = [];
   end
   
elseif method == 2

   % LIFT and BOX
   LIFT = zeros(d,S);
   BOX  = zeros(d,S);
   
   p = S;    % p = s^d;
   q = 1;
   
   for k = 1:d
      % LIFT
      LIFT(k,:) = mod(0:q:q*(S-1),S-1);
      
      % BOX
      p = p/s;
      q = q*s;
      kk = 1:q;
      tmp = repmat(mod(kk-1,s)+1,p,1);
      BOX(k,:) = tmp(:)-1;
   end
   
   LIFT = mod( LIFT , D );
   LIFT(:,end) = D-1;
   
   if 1
      % Not completely ready. It works, but the pattern for n points is static.
      x = orthSamp(n, 1, x_L, x_U, -1, 0);
      X = 1 + repmat((x-1)*D,1,S) + reshape(repmat(LIFT + BOX*M,n,1),d,N);
   else
      for k = 1:size(BOX,2)
         x = orthSamp(n, 1, x_L, x_U, -1, 0);
         X(:,(k-1)*n+1:k*n) = 1 + (x-1)*D + repmat( BOX(:,k)*M + LIFT(:,k) ,1,n);
      end
   end
   
end


if IntVals == 0
   X = repmat(x_L,1,N) + ((X-1) + rand(d,N)).*repmat(xD,1,N)/N;
end


% Possible values:
%
% d= 2. N = [1 4 8 9 12 16 18 20 24 25 27 28 32 36 40 44 45 48 49 50 52 54 56 60 63 64 68 72 75 76 80 81 84 88 90 92 96 98 99 100]
% d= 3. N = [1 8 16 24 27 32 40 48 54 56 64 72 80 81 88 96 104 108 112 120 125 128 135 136 144 152 160 162 168 176 184 189 192 200]
% d= 4. N = [1 16 32 48 64 80 81 96 112 128 144 160 162 176 192 208 224 240 243 256 272 288 304 320 324 336 352 368 384 400 405]
% d= 5. N = [1 32 64 96 128 160 192 224 243 256 288 320 352 384 416 448 480 486 512 544 576 608 640 672 704 729 736 768 800 832]
% d= 6. N = [1 64 128 192 256 320 384 448 512 576 640 704 729 768 832 896 960 1024 1088 1152 1216 1280 1344 1408 1458 1472 1536 1600]
% d= 7. N = [1 128 256 384 512 640 768 896 1024 1152 1280 1408 1536 1664 1792 1920 2048 2176 2187 2304 2432 2560 2688 2816 2944 3072 3200]
% d= 8. N = [1 256 512 768 1024 1280 1536 1792 2048 2304 2560 2816 3072 3328 3584 3840 4096 4352 4608 4864]
% d= 9. N = [1 512 1024 1536 2048 2560 3072 3584 4096 4608 5120 5632 6144 6656 7168 7680 8192 8704 9216 9728 10240]
%
% d=10. N = [1 1024 2048 3072 4096 5120 6144 7168 8192 9216 10240 11264 12288]
% d=11. N = [1 2048 4096 6144 8192 10240 12288 14336 16384 18432]
% d=12. N = [1 4096 8192 12288 16384]
% d=13. N = [1 8192 16384 24576]
% d=14. N = [1 16384 32768]
% d=15. N = [1 32768];



% MODIFICATION LOG:
%
% 080410 hkh Written 
% 080608 hkh Add PriLev, only print if PriLev > 1 or PriLev > 2
% 080608 hkh Remove check for nargin, demand all input
% 080714 nhq Revised most of the code, added new method.
% 080715 nhq Possible to return design as integer values.
% 090824 hkh Minor mlint revision
