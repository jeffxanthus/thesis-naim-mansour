% function ellipsePlot(ellipsesX, ellipsesY, ellipsesZ, axesA, axesB, axesC)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc., $Release: 6.1.0$
% Written April 10, 2008.    Last modified April 10, 2008.

function ellipsePlot(ellipsesX, ellipsesY, ellipsesZ, axesA, axesB, axesC)

[d,ellipseCount] = size(ellipsesX);

if d == 2
    if nargin == 1
        ellipsAxes = ones(d,1);
    elseif length(ellipsesY) == 1
        ellipsAxes = ellipsesY*ones(d,1);
    else
        ellipsAxes = ellipsesY(:);
    end
elseif d == 3
    if nargin < 6
        ellipsesAxes = [];
        axesC = [];
    end

    if nargin == 1
        if size(ellipsesX,1) == 4
            ellipsesAxes = ellipsesX(4,:);
        elseif size(ellipsesX,1) == 6
            axesA = ellipsesX(4,:);
            axesB = ellipsesX(5,:);
            axesC = ellipsesX(6,:);
        end
        ellipsesY = ellipsesX(2,:);
        ellipsesZ = ellipsesX(3,:);
        ellipsesX = ellipsesX(1,:);
    elseif nargin == 2
        s = size(ellipsesY);
        if all( s == 1 )
            ellipsesAxes = ellipsesY;
        elseif all( s == [1 3] )  ||  all( s == [3 1] )
            axesA = ones(1,ellipseCount)*ellipsesY(1);
            axesB = ones(1,ellipseCount)*ellipsesY(2);
            axesC = ones(1,ellipseCount)*ellipsesY(3);
        elseif all( s == [3 ellipseCount] )
            axesA = ellipsesY(1,:);
            axesB = ellipsesY(2,:);
            axesC = ellipsesY(3,:);
        end
        ellipsesY = ellipsesX(2,:);
        ellipsesZ = ellipsesX(3,:);
        ellipsesX = ellipsesX(1,:);
    elseif nargin == 3
        ellipsesAxes = 1;
    elseif nargin == 4
        s = size(axesA);
        if all( s == 1 )
            ellipsesAxes = axesA;
        elseif all( s == [1 3] )  ||  all( s == [3 1] )
            axesB = ones(1,ellipseCount)*axesA(2);
            axesC = ones(1,ellipseCount)*axesA(3);
            axesA = ones(1,ellipseCount)*axesA(1);
        elseif all( s == [ellipseCount 3] )  ||  all( s == [3 ellipseCount] )
            axesB = axesA(2,:);
            axesC = axesA(3,:);
            axesA = axesA(1,:);
        end
    end

    if isempty(axesC)
        if isempty(ellipsesAxes)
            axesA = ones(1,ellipseCount);
        elseif length(ellipsesAxes) == 1
            axesA = ellipsesAxes*ones(1,ellipseCount);
            %     elseif length(ellipsesAxes) == 3
            %         axesA = ellipsesAxes*ones(1,ellipseCount);
        elseif length(ellipsesAxes) == ellipseCount
            axesA = ellipsesAxes;
        end
        axesB = axesA;
        axesC = axesA;
    end
end

% set up basic plot
figure
axis equal
hold on

if d == 2;
    plot(ellipsesX(1,:),ellipsesX(2,:),'*');
    for i = 1:ellipseCount
        rectangle('Position',[ellipsesX(:,i)-ellipsAxes ; 2*ellipsAxes],'Curvature',[1 1]);
    end
    maxAxs = max(ellipsAxes,[],2);
    set(gca,'XLim',[floor(min(ellipsesX(1,:)))-maxAxs(1) ceil(max(ellipsesX(1,:)))+maxAxs(1)]);
    set(gca,'YLim',[floor(min(ellipsesX(2,:)))-maxAxs(1) ceil(max(ellipsesX(2,:)))+maxAxs(1)]);
elseif d == 3
    for i=1:ellipseCount
        ellipsoid(ellipsesX(i),ellipsesY(i),ellipsesZ(i),axesA(i),axesB(i),axesC(i))
    end
    colormap([0.7 0 0])
    view(3);
    rotate3d on;
    maxAxs = max([axesA ; axesB ; axesC],[],2);
    set(gca,'XLim',[floor(min(ellipsesX))-maxAxs(1) ceil(max(ellipsesX))+maxAxs(1)]);
    set(gca,'YLim',[floor(min(ellipsesY))-maxAxs(2) ceil(max(ellipsesY))+maxAxs(2)]);
    set(gca,'ZLim',[floor(min(ellipsesZ))-maxAxs(3) ceil(max(ellipsesZ))+maxAxs(3)]);
end

hold off

% x_L = floor(min(X,[],2));
% x_U = ceil(max(X,[],2));
% set(gca,'XLim',[x_L(1)-1 x_U(1)+1]);
% set(gca,'YLim',[x_L(2)-1 x_U(2)+1]);
% set(gca,'ZLim',[x_L(3)-1 x_U(3)+1]);


% MODIFICATION LOG:
%
% 080410 hkh Written 
