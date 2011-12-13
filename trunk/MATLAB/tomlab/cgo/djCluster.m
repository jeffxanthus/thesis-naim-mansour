% function [v,DistGrp] = djCluster(X,onB,epsX,PriLev,IntVars,Reals)
%
% IV    Logical vector of length Prob.N, 1 = integer variable
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2007 by Tomlab Optimization Inc., $Release: 4.0.0$
% Written Nov 24, 2006. Last modified March 6, 2007.

function [v,DistGrp] = djCluster(X,onB,epsX,PriLev,IntVars,Reals)

if nargin < 6
   Reals = [];
   if nargin < 5
      IntVars = [];
      if nargin < 4
         PriLev = [];
         if nargin < 3
            epsX = [];
            if nargin < 2
               onB = [];
               if nargin < 1
                  error('djCluster needs one parameter, the matrix X');
               end
            end
         end
      end
   end
end
if isempty(PriLev), PriLev = 1; end
if isempty(epsX),   PriLev = 1E-4; end

[d,m]    = size(X);

if isempty(Reals)
   IV          = ones(d,1);
   IV(IntVars) = 0;
   Reals       = find(IV);
end

if isempty(onB), onB = zeros(m,1); end

v        = ones(1,m);
dsqrt    = sqrt(length(Reals));
delta    = zeros(m,2);
DistGrp  = zeros(m,1);
if isempty(IntVars)
   for i = 1:m-1
       delta(i,1) = tomsol(30,X(:,i),X(:,i+1))/dsqrt;
   end
else
   for i = 1:m-1
       delta(i,1) = tomsol(30,X(Reals,i),X(Reals,i+1))/dsqrt;
       delta(i,2) = sum(X(IntVars,i)~=X(IntVars,i+1));
   end
end
Group = 1;
LastG = 1;
if PriLev > 0
   fprintf('Target Group Delta(i) ');
   if ~isempty(IntVars)
      fprintf(' IV ');
   end
   fprintf(' Criterion GrpDistance    x(1),...,x(d)\n');
   fprintf('%4d %5d  %9.5f ', 1,1,delta(1,1));
   if ~isempty(IntVars)
      fprintf('%3d ',delta(1,2));
   end
   fprintf('%10.5f %9.5f ', NaN,0);
   xprint(X(:,1),'x:')
end
for i = 2:m
    if delta(i-1,2) >= 1 
       %HKH special - Integer values have changed
       C = 999;
    elseif delta(i,1) > 0.1 & delta(i-1,1) > 0.1
       C = 100;
    elseif delta(i,1) > 0.0005
       C = delta(i-1,1)/delta(i,1);
    elseif i>=3 & delta(i-1,1) > 0.0005
       C = delta(i-1,1)/max(delta(i-2,1),0.0005);
    elseif i==2 & delta(1,1) > 0.1 & delta(2,1) < 0.0005
       C = 100;
    else
       C = 0;
    end
    if i==m & C < 12 & delta(m-1,1) > 0.1
       %HKH special
       C = 100+C;
    end
    if C < 12 & delta(i-1,1) >= 0.1 & onB(i) ~= onB(i-1)
       %HKH special, check for 10% change in component on bound 
       for k = 1:length(Reals)
           j = Reals(k);
           if X(j,i)   >= 1-epsX & X(j,i-1) <= 0.9 | ...
              X(j,i)   <= epsX   & X(j,i-1) >= 0.1 | ...
              X(j,i-1) >= 1-epsX & X(j,i)   <= 0.9 | ...
              X(j,i-1) <= epsX   & X(j,i)   >= 0.1 
              C = 200+C;
           end
       end
    end
    if C >= 12
       Group = Group + 1;
       DistGrp(i) = 0; 
       LastG      = i;
    else
       DistGrp(i) = tomsol(30,X(Reals,i),X(Reals,LastG))/dsqrt;
    end
    v(i) = Group;
    if PriLev > 0
       fprintf('%4d %5d  %9.5f ', i,Group,delta(i,1));
       if ~isempty(IntVars)
          fprintf('%3d ',delta(i,2));
       end
       fprintf('%10.5f %9.5f ',C,DistGrp(i));
       %fprintf('%4d %5d  %9.5f %10.5f %9.5f ', ...
       %        i,Group,delta(i), C,DistGrp(i));
       xprint(X(:,i),'x:')
    end
end


% MODIFICATION LOG
%
% 061219 hkh  Add input epsX. Add C to 100 and 200, to see the original C value
% 070306 hkh  Special treatment of mixed-integer problems
