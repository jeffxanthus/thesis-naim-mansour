%
% function [X,nCon] = randomDesigns(method,x_L,x_U,percent,nSample,dist,X,...
%                                   nTrial,Prob,nCon,SCALE,PriLev)
% This function sample feasible points inside box-bounds. To avoid sampling
% points very close, and to spread out the points over the feasible space,
% a surrounding to each point is considered infeasible.
%
% Possible methods:  1. Sphere,  2. Ellipse,  3. Rectangle
%
% It is possible to sample only feasible points with respect to linear and
% nonlinear constraints. This requires a Prob structure as input argument.
% Instead of method = 1, 2 or 3, set method to -1, -2 or -3 respectively.
%
% NOTE: The calculations to see if ellipsoids are overlapping is much more
%       demanding than checking for spheres overlap. Therefore, if the box
%       formed by x_L and x_U is near-quadratic, one should use the Sphere
%       method instead since the results will be very similair.
%
%       The main point with Ellipsoids is to account for spaces with a big
%       difference in range for the dimensions.
%
%       If SCALE is set, all ranges are 1. This makes the Ellipsoid Method 
%       useless, as it will result in the same sample as Circle Method.
%
%       More general, if all ranges are equal (not necessary 1), the same
%       argumentation is valid. Therefore, if Ellipsoid is chosen but all
%       the ranges are equal, automatically switch to Circle instead.
%
% function [X,nCon] = randomDesigns(method,x_L,x_U,percent,nSample,dist,X,...
%                                   nTrial,Prob,nCon,SCALE,PriLev)
%
% Input parameters:
% method        Choose between   1 "Sphere", 2 "Ellipse" or 3 "Rectangle".
% x_L, x_U      Lower and Upper bounds for the variables x.
% percent       Percentage of the variable space to surround each point.
% nSample       Number of points to be sampled.
% dist          Distance, used to override default surrounding.
% X             Matrix of already sampled points.
% nTrial        maximum number of trial points randomly generated
%               for each new point to sample.
%
%           ONLY NEEDED FOR THE CONSTRAINED VERSIONS
% Prob          Problem structure
% nCon          Number of constraint evaluations.
% SCALE         0 - Original search space
%               1 - Transform search space to unit cube
%
% Output parameters:
% X             Matrix of sampled points, including eventual input X.
% nCon          Number of constraint evaluations.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written April 10, 2008.    Last modified Aug 13, 2009.

function [X,nCon] = randomDesigns(method,x_L,x_U,percent,nSample,dist,X,...
                                  nTrial,Prob,nCon,SCALE,PriLev)

if nargin < 12
   PriLev = 2;
   if nargin < 11
      SCALE = 0;
      if nargin < 10
         nCon = 0;
         if nargin < 9
            Prob = [];
            if nargin < 8
               nTrial = [];
               if nargin < 7
                  X = NaN;
                  if nargin < 6
                     dist = [];
                     if nargin < 5
                        nSample = 2*length(x_L);
                     end
                  end
               end
            end
         end
      end
   end
end

% Number of random testpoints each iteration
if isempty(nTrial)
    MAXTEST  = 1000;
else
    MAXTEST = nTrial;
end

INTMethod = 2;    % Treat IntVars in different ways.
                  % 1 - Deals with IntVars after generation of points.
                  % 2 - Round IntVars BEFORE distance check to sampled points.
PLOTFLAG  = 0;    % If set to 1, a plot is shown (only if d = 2 or 3)
BOXSIZE   = 2;    % Factor to even out the percent option. If not set,
                  % the percent variable should be set twice as high
                  % to get a similair surrounding compared to Spheres.

d  = length(x_L);
xD = x_U - x_L;

Constrained = 0;
if method < 0
    method = abs(method);
    if isempty(Prob)
       if PriLev > 1
          fprintf('A Prob structure is needed to find feasible points.\n')
          fprintf('Method switched to unconstrained version.\n')
       end
    else
        Constrained = 1;
    end
end

if all(xD(1) == xD(2:end))
   % Unit cube space, unnecessary to use Ellipsoid Method.
   if method == 2
      if PriLev > 1
         fprintf(' Range xD equal in all dimensions, switch')
         fprintf(' method from Ellipsoid to Circle.\n')
      end
      method = 1;
   end
end

if ~isempty(Prob)
   IntVars = Prob.MIP.IntVars;
else
   IntVars = [];
end


% N looping parameter
N = nSample;
if isnan(X)
    X  = [x_L + xD/2,zeros(d,nSample-1)];             % center point
    k  = 1;
    %NHQ Test this point for feasibility as well???
elseif isempty(X)
    xStart = x_L + rand(d,1).*xD;                     % random point
    % Round any IntVars meanwhile generating points
    if INTMethod == 2  &&  ~isempty(IntVars)
       xStart(IntVars,:) = round(xStart(IntVars,:));
    end
    if Constrained
        for i = 1:MAXTEST
            [flag,nCon] = isFeasible(xStart,Prob,SCALE,nCon);
            if flag == 1    % X feasible, quit loop.
                break
            else            % Generate new random point
                xStart = x_L + rand(d,1).*xD;
                if INTMethod == 2  &&  ~isempty(IntVars)
                   xStart(IntVars,:) = round(xStart(IntVars,:));
                end
            end
        end
        if i >= MAXTEST  &&  flag == 0
           if PriLev > 1
              fprintf(' Couldn''t find any feasible point.\n')
              fprintf(' Continue with last tested infeasible X\n')
           end
        end
    end
    X  = [xStart,zeros(d,nSample-1)];
    k  = 1;
else
    k  = size(X,2);           % Number of given points
    X(d,k+1:k+nSample) = 0;   % Initialize
    N  = N+1;                 % Increase N, because looping to nSample-1
end

% If dist == [], use xD as starting point.
% For Sphere, the diameter is set to norm(xD)
diam  = dist;
if isempty(dist)
    % Use xD as range in each dimension respectively
    dist = xD;
    %diam = norm(xD);
    diam = 2*min(xD);
elseif isinf(dist)
    % Use norm(xD) as range in all dimensions (square)
    dist = norm(xD)*ones(d,1);
    diam = 2*norm(xD);
elseif length(dist) ~= d
    % Length mismatch, use dist(1) in all dimensions (square)
    dist = dist(1)*ones(d,1);
    diam = 2*dist(1);
else
    dist = dist(:);
end

if method == 1;
    % Put a sphere with radius 'minDist' around each sampled point
    minDist = percent/100*diam;
elseif method == 2
    % Put an ellipse with axes 'ellipsAxes' around each sampled point
    ellipsAxes = percent/100*dist;
else    %elseif method == 3
    % Put a rectangle with sides 'box' around each sampled point
    box = BOXSIZE*percent/100*dist;
end

% Find n points not overlapping, or as spread out as possible
for i=1:N-1
    j = 0;
    mInfeasible = 0;
    while j < MAXTEST
        j    = j+1;
        xNew = x_L + rand(d,1).*xD;
        % Round any IntVars meanwhile generating points
        if INTMethod == 2  &&  ~isempty(IntVars)
           xNew(IntVars,:) = round(xNew(IntVars,:));
        end
        % Accept only feasible points??
        if Constrained
            [flag,nCon] = isFeasible(xNew,Prob,SCALE,nCon);
            if flag == 0
                % xNew not feasible, skip rest of loop.
                continue
            end
        end

        if method == 1
            % Check if Spheres overlap
            MD   = min(tomsol(30,xNew,X(:,1:k+i-1)));
            % test if point is OK
            if MD >= minDist
                break
            elseif MD > mInfeasible
                mInfeasible = MD;
                xBest = xNew;
            end
        else
            if method == 2
                % Check if Ellipsoids overlap
                MD = checkEllipseOverlap(xNew,X(:,1:k+i-1),ellipsAxes);
            else    %elseif method == 3
                % Check if Rectangles overlap
                MD = checkBoxes(xNew,X(:,1:k+i-1),box);
            end
            % test if point is OK
            if MD == 0
                break
            elseif MD > mInfeasible
                mInfeasible = MD;
                xBest = xNew;
            end
        end
    end
    if j >= MAXTEST
        xNew = xBest;   % No "feasible" point found, add best found instead
    end
    X(:,k+i) = xNew;
end


if ~isempty(IntVars)  &&  Constrained == 0  &&  INTMethod == 1
   % Use experimental design values rounded to nearest integer point
   X  = daceInts(X,xD,x_U,IntVars);
   if size(X,2) ~= k+nSample
      % fprintf('Removed %d duplicate points when sampling Integer points.\n',nSample-size(X,2))
   end
end


% Plot the sampled points if PLOTFLAG is set
if ( d == 2  ||  d == 3 )  &&  PLOTFLAG == 1
   if method == 1
      ellipsePlot(X,0.5*minDist)
   elseif method == 2
      ellipsePlot(X,ellipsAxes)
   else
      rectanglePlot(X,box)
   end
end


% ====================================================================
function MD = checkEllipseOverlap(xNew,X,ellipsAxes)
% ====================================================================
% Check if ellipsAxes around xNew will overlap any
% of the ellipsAxeses around points in X.
k  = size(X,2);
d2 = (repmat(xNew,1,k) - X).^2;
n  = sqrt(sum(d2));         % n = tomsol(30,xNew,X);
r  = sqrt(sum(d2.*repmat(ellipsAxes.^2,1,k)))./n;

% d  = repmat(xNew,1,k) - X;
% n  = sqrt(sum(d.^2));       % n = tomsol(30,xNew,X);
% r  = sqrt(sum((d.*repmat(ellipsAxes,1,k)).^2))./n;

if any( n == 0 )
    MD = inf;
    return
end

flag = 2*r <= n;

if all(flag)
    MD = 0;
else
    MD = min(n);
end


% ====================================================================
function MD = checkBoxes(xNew,X,box)
% ====================================================================
% Check if box around xNew will overlap
% any of the boxes around points in X.

if any(ismember(X',xNew','rows'))
    MD  = inf;
    return
end

k = size(X,2);
X0 = X - repmat(0.5*box,1,k);
X1 = X0 + repmat(box,1,k);
x0 = xNew - 0.5*box;
x1 = x0 + box;

flag = repmat([x1 ; -x0],1,k) < [X0 ; -X1];
if all(any(flag))
    MD = 0;
else
    MD = min(tomsol(30,xNew,X));
end


% ========================================================================
function [flag,nCon] = isFeasible(xNew,Prob,SCALE,nCon)
% ========================================================================
% Check if points in xNew are feasible with respect to given constrints

dLin = length(Prob.b_L);
dCon = length(Prob.c_L);
bTol = Prob.optParam.bTol;      % Linear constraint feasibility tolerance
cTol = Prob.optParam.cTol;      % Constraint feasibility tolerance

%Scale back to original space before checking constraints
if SCALE
    oNew = tomsol(9, Prob.x_L, xNew, Prob.x_U-Prob.x_L);
else
    oNew = xNew;
end

[d,n] = size(xNew);
flag  = zeros(n,1);

% loop and check all points
for i = 1:n
    cErr = 0;
    if dLin > 0
        L = Prob.A*oNew(:,i);
%        keyboard
        cErr = cErr+sum(max(0,max(Prob.b_L-bTol-L,L-bTol-Prob.b_U)));
    end
    if dCon > 0 && cErr == 0
        %C = nlp_c(oNew(:,i), Prob, varargin{:});
        C = nlp_c(oNew(:,i), Prob);
        nCon = nCon + 1;
        cErr = cErr+sum(max(0,max(Prob.c_L-cTol-C,C-cTol-Prob.c_U)));
    end
    if cErr == 0
        flag(i) = 1;
    end
end

% MODIFICATION LOG:
%
% 080410 nhq Written 
% 080410 hkh Revised, was awful mess, added Prilev 
% 080709 nhq Methods now handle both IntVars and Constraints
% 080916 nhq Corrected bug in function "checkEllipseOverlap"
% 090813 med repmat call corrected
% 090918 hkh Change call to daceInts, now external
