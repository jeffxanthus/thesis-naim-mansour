% addCloseX.m
%
% function xLoc = addCloseX(xLoc,X,x_L,x_U,IntVars)
%
% addCloseX finds mixed-integer x value as close as possible to xLoc
%
% INPUT:
%
% M       Number of trial points
%
% xLoc    Local minimum, already part of X
% X       Set X of sampled points
% x_L     Lower bounds for each element in x.
% x_U     Upper bounds for each element in x.
% IntVars Indices for integer variables in x
%
%
% OUTPUT:
% xLoc    New x, closest to xLoc, but not in X

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written June 4, 2008.   Last modified August 24, 2009.

function xLoc = addCloseX(xLoc,X,x_L,x_U,IntVars)

% Try to find a neighbour integer point
d     = length(x_L);
L     = length(IntVars);
m     = 0;
if L == d
   x = xLoc;
   % Check size 1 changes of xLoc
   for i = 1:d
       x(i) = xLoc(i)+1;
       m    = m+1;
       if x(i) <= x_U(i)
          dXsny    = min(tomsol(30,x,X));
          PrintdXsny(i,0,0,x,dXsny)
          if dXsny > 0, 
             xLoc = x;
             return; 
          end
       end
       x(i) = xLoc(i)-1;
       m    = m+1;
       if x(i) >= x_L(i)
          dXsny    = min(tomsol(30,x,X));
          PrintdXsny(-i,0,0,x,dXsny)
          if dXsny > 0, 
             xLoc = x;
             return; 
          end
       end
       x(i) = xLoc(i);
   end
   % Check size 2 changes of xLoc
   for i = 1:d-1
       x(i) = xLoc(i)+1;
       if x(i) <= x_U(i)
          for j = i+1:d
              x(j) = xLoc(j)+1;
              m    = m+1;
              if x(j) <= x_U(j)
                 dXsny    = min(tomsol(30,x,X));
                 PrintdXsny(i,j,0,x,dXsny)
                 if dXsny > 0, 
                    xLoc = x;
                    return; 
                 end
              end
              x(j) = xLoc(j)-1;
              m    = m+1;
              if x(j) >= x_L(j)
                 dXsny    = min(tomsol(30,x,X));
                 PrintdXsny(i,-j,0,x,dXsny)
                 if dXsny > 0, 
                    xLoc = x;
                    return; 
                 end
              end
              % Reset back to xLoc(j)
              x(j) = xLoc(j);
          end
       end
       x(i) = xLoc(i)-1;
       if x(i) >= x_L(i)
          for j = i+1:d
              x(j) = xLoc(j)+1;
              m    = m+1;
              if x(j) <= x_U(j)
                 dXsny    = min(tomsol(30,x,X));
                 PrintdXsny(-i,j,0,x,dXsny)
                 if dXsny > 0, 
                    xLoc = x;
                    return; 
                 end
              end
              x(j) = xLoc(j)-1;
              m    = m+1;
              if x(j) >= x_L(j)
                 dXsny    = min(tomsol(30,x,X));
                 PrintdXsny(-i,-j,0,x,dXsny)
                 if dXsny > 0, 
                    xLoc = x;
                    return; 
                 end
              end
              % Reset back to xLoc(j)
              x(j) = xLoc(j);
          end
       end
       % Reset back to xLoc(i)
       x(i) = xLoc(i);
   end
   % Check size 3 changes of xLoc
   for i = 1:d-2
       x(i) = xLoc(i)+1;
       if x(i) <= x_U(i)
          for j = i+1:d-1
              x(j) = xLoc(j)+1;
              if x(j) <= x_U(j)
                 for k = j+1:d
                     x(k) = xLoc(k)+1;
                     if x(k) <= x_U(k)
                        m    = m+1;
                        dXsny    = min(tomsol(30,x,X));
                        PrintdXsny(i,j,k,x,dXsny)
                        if dXsny > 0, 
                           xLoc = x;
                           return; 
                        end
                     end
                     x(k) = xLoc(k)-1;
                     if x(k) >= x_L(k)
                        m    = m+1;
                        dXsny    = min(tomsol(30,x,X));
                        PrintdXsny(i,j,-k,x,dXsny)
                        if dXsny > 0, 
                           xLoc = x;
                           return; 
                        end
                     end
                     % Reset back to xLoc(k)
                     x(k) = xLoc(k);
                 end
              end
              x(j) = xLoc(j)-1;
              if x(j) >= x_L(j)
                 for k = j+1:d
                     x(k) = xLoc(k)+1;
                     if x(k) <= x_U(k)
                        m    = m+1;
                        dXsny    = min(tomsol(30,x,X));
                        PrintdXsny(i,-j,k,x,dXsny)
                        if dXsny > 0, 
                           xLoc = x;
                           return; 
                        end
                     end
                     x(k) = xLoc(k)-1;
                     if x(k) >= x_L(k)
                        m    = m+1;
                        dXsny    = min(tomsol(30,x,X));
                        PrintdXsny(i,-j,-k,x,dXsny)
                        if dXsny > 0, 
                           xLoc = x;
                           return; 
                        end
                     end
                     % Reset back to xLoc(k)
                     x(k) = xLoc(k);
                 end
              end
              % Reset back to xLoc(j)
              x(j) = xLoc(j);
          end
       end
       x(i) = xLoc(i)-1;
       if x(i) >= x_L(i)
          for j = i+1:d-1
              x(j) = xLoc(j)+1;
              if x(j) <= x_U(j)
                 for k = j+1:d
                     x(k) = xLoc(k)+1;
                     if x(k) <= x_U(k)
                        m    = m+1;
                        dXsny    = min(tomsol(30,x,X));
                        PrintdXsny(-i,j,k,x,dXsny)
                        if dXsny > 0, 
                           xLoc = x;
                           return; 
                        end
                     end
                     x(k) = xLoc(k)-1;
                     if x(k) >= x_L(k)
                        m    = m+1;
                        dXsny    = min(tomsol(30,x,X));
                        PrintdXsny(-i,j,-k,x,dXsny)
                        if dXsny > 0, 
                           xLoc = x;
                           return; 
                        end
                     end
                     % Reset back to xLoc(k)
                     x(k) = xLoc(k);
                 end
              end
              x(j) = xLoc(j)-1;
              if x(j) >= x_L(j)
                 for k = j+1:d
                     x(k) = xLoc(k)+1;
                     if x(k) <= x_U(k)
                        m    = m+1;
                        dXsny    = min(tomsol(30,x,X));
                        PrintdXsny(-i,-j,k,x,dXsny)
                        if dXsny > 0, 
                           xLoc = x;
                           return; 
                        end
                     end
                     x(k) = xLoc(k)-1;
                     if x(k) >= x_L(k)
                        m    = m+1;
                        dXsny    = min(tomsol(30,x,X));
                        PrintdXsny(-i,-j,-k,x,dXsny)
                        if dXsny > 0, 
                           xLoc = x;
                           return; 
                        end
                     end
                     % Reset back to xLoc(k)
                     x(k) = xLoc(k);
                 end
              end
              % Reset back to xLoc(j)
              x(j) = xLoc(j);
          end
       end
       % Reset back to xLoc(i)
       x(i) = xLoc(i);
   end
   fprintf('addCloseX: Failed to find ');
   fprintf('neighbour integer point, m== %d tries\n',m);
end

function PrintdXsny(d1,d2,d3,x,dXsny)
fprintf('d1 = %3d.',d1);
if d2 ~= 0 
   fprintf(' d2 = %3d.',d2);
end
if d3 ~= 0 
   fprintf(' d3 = %3d.',d3);
end

fprintf('dXsny = %f:',dXsny);
xprinti(x)

% MODIFICATION LOG:
%
% 080604  hkh  Written
% 080615  hkh  New algorithm, find 1,2 or 3-step nearest point
