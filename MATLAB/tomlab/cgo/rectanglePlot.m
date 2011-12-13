% function rectanglePlot(X,box)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc., $Release: 6.1.0$
% Written April 10, 2008.    Last modified April 10, 2008.

function rectanglePlot(X,box)

[d,nVars] = size(X);

if nargin == 1
    box = ones(d,1);
elseif length(box) == 1
    box = ones(d,1);
else
    box = box(:);
end

figure;
%clf
axis equal;
hold on

if d == 2;
    plot(X(1,:),X(2,:),'*');
    for i = 1:nVars
        rectangle('Position',[X(:,i) - 0.5*box ; box],'Curvature',[0 0]);
    end
elseif d == 3
    for i = 1:nVars
        X0 = X(:,i)-0.5*box;
        X1 = X0 + box;

        x1 = X0(1);   x2 = X1(1);
        y1 = X0(2);   y2 = X1(2);
        z1 = X0(3);   z2 = X1(3);

        Vertices = [x1, y1, z1;  % 1
            x2, y1, z1;  % 2
            x1, y1, z2;  % 3
            x2, y1, z2;  % 4
            x1, y2, z1;  % 5
            x2, y2, z1;  % 6
            x1, y2, z2;  % 7
            x2, y2, z2]; % 8

        Faces = [1, 2, 4, 3;  % back face
                 1, 2, 6, 5;  % bootom face
                 1, 3, 7, 5;  % left face
                 7, 8, 6, 5;  % front face
                 2, 4, 8, 6;  % right face
                 3, 4, 8, 7]; % upper face
        
        patch('Vertices',Vertices,'Faces', Faces, ...
              'FaceVertexCData',hsv(1),'FaceColor','flat');
    end
    colormap([0.7 0 0])
    view(3);
    rotate3d on;
end
box = 0.5*box;
set(gca,'XLim',[floor(min(X(1,:)))-box(1) ceil(max(X(1,:)))+box(1)]);
set(gca,'YLim',[floor(min(X(2,:)))-box(2) ceil(max(X(2,:)))+box(2)]);
if d == 3
    set(gca,'ZLim',[floor(min(X(3,:)))-box(3) ceil(max(X(3,:)))+box(3)]);
end

% MODIFICATION LOG:
%
% 080410 hkh Written 
