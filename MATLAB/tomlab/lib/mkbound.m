% Make bounds of same size and put default +-Inf, if missing
%
% function [Prob] = mkbound(Prob);
%
% INPUT:
%  Prob     Structure. Using:
%           Prob.x_0,Prob.A
%           Prob.x_L,Prob.x_U
%           Prob.b_L,Prob.b_U
%           Prob.c_L,Prob.c_U
%
% OUTPUT:
%  Prob     Structure. Changing:
%           Prob.x_L, Prob.x_U, Prob.b_L, Prob.b_U, Prob.c_L, Prob.c_U

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Mar 2, 1998.   Last modified Aug 13, 2009.

function [Prob] = mkbound(Prob)

n=max([length(Prob.x_0),length(Prob.x_L),length(Prob.x_U)]);

if n == 0, n=size(Prob.A,2); end

if isempty(Prob.x_L)
   Prob.x_L=-Inf*ones(n,1);
elseif length(Prob.x_L) < n
   Prob.x_L=[Prob.x_L(:); -Inf*ones(n-length(Prob.x_L),1)];
else
   Prob.x_L=Prob.x_L(:);
end

if isempty(Prob.x_U)
   Prob.x_U=Inf*ones(n,1);
elseif length(Prob.x_U) < n
   Prob.x_U=[Prob.x_U(:); Inf*ones(n-length(Prob.x_U),1)];
else
   Prob.x_U=Prob.x_U(:);
end

if isfield(Prob,'A')
   if ~isempty(Prob.A)
      m1=size(Prob.A,1);
   
      if size(Prob.A,2) ~= n
         fprintf('Size A %d %d. n = %d\n',size(Prob.A),n);
         error('mkbound: illegal size for A');
      end

      if isempty(Prob.b_L)
         Prob.b_L=-Inf*ones(m1,1);
      elseif length(Prob.b_L) ~= m1
         error('mkbound: illegal size on b_L');
      else
         Prob.b_L=Prob.b_L(:);
      end

      if isempty(Prob.b_U)
         Prob.b_U=Inf*ones(m1,1);
      elseif length(Prob.b_U) ~= m1
         error('mkbound: illegal size on b_U');
      else
         Prob.b_U=Prob.b_U(:);
      end
   end
end

m2=max(length(Prob.c_L),length(Prob.c_U));

if m2 > 0
   if isempty(Prob.c_L)
      Prob.c_L=-Inf*ones(m2,1);
   elseif length(Prob.c_L)~=m2
      Prob.c_L
      error('mkbound: illegal size on c_L');
   else
      Prob.c_L=Prob.c_L(:);
   end
   if isempty(Prob.c_U)
      Prob.c_U=Inf*ones(m2,1);
   elseif length(Prob.c_U)~=m2
      error('mkbound: illegal size on c_U');
   else
      Prob.c_U=Prob.c_U(:);
   end
end

% Generate a start vector x0 if empty, (and not LP problems)
if isempty(Prob.x_0) & Prob.probType~=8
   Prob.x_0=max(Prob.x_L,min(zeros(n,1),Prob.x_U));
   Prob.N=n;    % Set to correct length
end

% MODIFICATION LOG:
%
% 980930  hkh  Terrible sign bug: Setting upper bound x_U to -Inf
% 981106  hkh  -Inf on x_U on line 35.
% 981107  hkh  Change definition of n, use max of x_L,x_U and x_0
%              Generate a start vector x0 if empty. New Prob.N value
% 981111  hkh  Do not generate x_0 for LP problems
% 090813  med  mlint check