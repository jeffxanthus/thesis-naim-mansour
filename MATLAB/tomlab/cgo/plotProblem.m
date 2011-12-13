% function fig = plotProblem(Prob,X,F,figIN,MaxFLAG,u,v)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc., $Release: 6.1.0$
% Written April 10, 2008.    Last modified April 10, 2008.

function fig = plotProblem(Prob,X,F,figIN,MaxFLAG,u,v)


d = Prob.N;

if ~ismember(d,[1 2])
   fig = [];
   return
end

bigM = 1000;

if d == 1
   gridSize = 1000;
else
   gridSize = 151;
end

% Info from Prob
prob_type   = Prob.probFile;
prob_number = Prob.P;
prob_name   = Prob.Name;
x_opt       = Prob.x_opt;
f_opt       = Prob.f_opt;
Func_f      = Prob.FUNCS.f;

if nargin == 1
   X = [];
   F = [];
elseif nargin == 2
   F = [];
end
if nargin <= 3
   figIN = [];
end
if nargin <= 4
   MaxFLAG = 0;
end
if nargin <= 5
   x_L     = Prob.x_L;
   x_U     = Prob.x_U;
   
   if isempty(x_L)
      x_L = -inf(d,1);
   end
   if isempty(x_U)
      x_U = inf(d,1);
   end
   
   idxL    = any([isinf(x_L)' ; isnan(x_L)'])';
   idxU    = any([isinf(x_U)' ; isnan(x_U)'])';
   x_L(idxL) = max(-bigM,Prob.x_min(idxL)');
   x_U(idxU) = min( bigM,Prob.x_max(idxU)');
%    x_L(idxL) = -bigM;
%    x_U(idxU) =  bigM;
%    x_L = max(x_L,Prob.x_min(:));
%    x_U = min(x_U,Prob.x_max(:));
   
   u       = linspace(x_L(1),x_U(1),gridSize);
   if d == 2
       v  = linspace(x_L(2),x_U(2),gridSize);
   else
       v = 1;
   end
end

w = zeros(length(v),length(u));

%keyboard

for i = 1:length(u)
   if d == 1
      x = u(i);
      w(i) = eval([Func_f '(x,Prob)']);
      continue
   end
   for j = 1:length(v)
      x = [u(i) v(j)];
      ff = eval([Func_f '(x,Prob)']);
      if ~isreal(ff)
          ff = inf;
      end
      w(j,i) = ff;
   end
end

I = w > 1E4;
if any(any(I))  &&  ~all(all(I))
   maxw = max(max(w(~I)));
   w(I) = maxw + 30;
end

if ~isempty(X)  &&  isempty(F)
   for k = 1:size(X,2)
      F(k) = eval([Func_f '(X(:,k),Prob)']);
   end
end
I = F > 1E4;
if any(I)
   maxw = max(max(w));
   F(I) = maxw + 30;
end
MAXW = max(max(w));

if isempty(figIN)
   fig = figure();
   set(fig,'Name',[prob_type ' ' int2str(prob_number) ' - ' prob_name]);
else
   figure(figIN)
end

if MaxFLAG
   w     = -w;
   F     = -F;
   f_opt = -f_opt;
end

if d == 1
    plot(u,w)
    hold;
else
    mesh(u,v,w)
    rotate3d on
    hold;
end

for k = 1:size(x_opt,1)
   if d == 1
      plot(x_opt(k),f_opt,'r*','MarkerSize',18)
   else
      plot3(x_opt(k,1),x_opt(k,2),f_opt,'r*','LineWidth',2,'MarkerSize',18)
      % plot box around x_opt. Size of box = 0.1*x_D
      plotBox(Prob,x_opt(k,:),x_L,x_U,0.1,MaxFLAG);
   end
end

if ~isempty(X);
   if d == 1
      plot(X,F,'k.','MarkerSize',15)
   else
      plot3(X(1,:),X(2,:),F,'k.','MarkerSize',15)
   end
end

% xlim([x_L(1),x_U(1)])
% ylim([x_L(2),x_U(2)])
plotConstraints(Prob,u,v,MAXW,MaxFLAG)
plotLinearConstraints(Prob,MAXW,MaxFLAG,x_L,x_U);
%keyboard
xlim([x_L(1),x_U(1)])
if d == 2
   ylim([x_L(2),x_U(2)])
end


function plotBox(Prob,x_opt,x_L,x_U,frac,MaxFLAG)
% plot box around x_opt. Size of box = 0.1*x_D

smallGrid = 25;

Func = Prob.FUNCS.f;
xBox = frac*(x_U-x_L);

xB_L = x_opt(:) - 0.5*xBox;
xB_U = x_opt(:) + 0.5*xBox;
idxL = xB_L > x_L;
idxU = xB_U < x_U;
xB_L = max(xB_L,x_L);
xB_U = min(xB_U,x_U);

xB   = linspace(xB_L(1),xB_U(1),smallGrid);
yB   = linspace(xB_L(2),xB_U(2),smallGrid);
BF   = zeros(4,smallGrid);

for k = 1:smallGrid
    BF(1,k) = eval([Func '([xB(k) xB_L(2)],Prob)']);
    BF(2,k) = eval([Func '([xB(k) xB_U(2)],Prob)']);
    BF(3,k) = eval([Func '([xB_L(1) yB(k)],Prob)']);
    BF(4,k) = eval([Func '([xB_U(1) yB(k)],Prob)']);
end
if MaxFLAG
   BF     = -BF;
end
if idxU(2)
   plot3(xB,xB_U(2)*ones(25,1),BF(2,:),'r','LineWidth',2)
end
if idxL(2)
   plot3(xB,xB_L(2)*ones(25,1),BF(1,:),'r','LineWidth',2)
end
if idxU(1)
   plot3(xB_U(1)*ones(25,1),yB,BF(4,:),'r','LineWidth',2)
end
if idxL(1)
   plot3(xB_L(1)*ones(25,1),yB,BF(3,:),'r','LineWidth',2)
end



function plotConstraints(Prob,u,v,MAXW,MaxFLAG)

probType = 'none';
if isequal(Prob.probFile,'glc_prob')
    probType = 'glc_prob';
elseif isequal(Prob.probFile,'glcIP_prob')
    probType = 'glcIP_prob';
elseif isequal(Prob.probFile,'minlp_prob') && Prob.P == 4
    probType = 'glcIP_prob';
elseif isequal(Prob.probFile,'con_prob')
    probType = 'con_prob';
end

if ~isequal(probType,'none')
   if isfield(Prob,'SCALE')
      if Prob.SCALE == 1
         u  = Prob.xL(1) + u*Prob.xD(1);
         v  = Prob.xL(2) + v*Prob.xD(2);
      end
   end

   cx = [];
   
   switch probType
       
   case {'glc_prob'}
       
       if Prob.P == 1  ||  Prob.P == 2
           %         cx = asin(2*(sin(2*pi*v)).^2)/(4*pi);
           %         cy = v;
           x = 0:0.01:0.25;
           y = sqrt(1/64-(x-0.125).^2);
           cx = zeros(2*4*5,length(x));
           cy = zeros(2*4*5,length(x));
           i = 0;
           for k = -2:1
               for l = -2:2
                   i = i+2;
                   cx([i-1 i],:) = [x ; x]+k*0.5;
                   cy([i-1 i],:) = [y ;-y]+l*0.5;
               end
           end
       elseif Prob.P == 3
           cx = [u ; u ; 0.2*(v-50).^2+55];
           cy = [700./u ; 5*(u/25).^2 ; v];
           %      cy = [700./u ; (u/25).^2 ; v];
       elseif Prob.P == 7
           cx = (0:0.05:1);
           cy = sqrt(1-cx.^2);
           cx = [cx ; cx];
           cy = [cy ; -cy];
       elseif Prob.P == 8
           cx = [u ; u ];
           cy = [700./u ; 5*(u/25).^2 ];
       elseif Prob.P == 9
           cx = [u ; u ; 0.2*(v-50).^2+55];
           cy = [700./u ; 5*(u/25).^2 ; v];
       elseif Prob.P == 10
           cx = u;
           cy = 700./u;
       elseif Prob.P == 11
           cx = 8.62*v.^3;
           cy = v;
       elseif Prob.P == 12
           %OBS!! Linear approximation of very tricky constraint.
           %      Really good approx. though for visualization.
           k1 = -1.735294117647059;
           k2 =  1.720338983050848;
           m1 =  1.604117647058823;
           m2 = -1.540508474576272;
           cx = u;
           cy = max(k1.*u+m1,k2.*u+m2);
       elseif Prob.P == 17  % TP 9
           cx = [u ; u];
           cy = [2*u.^4-8*u.^3+8*u.^2+2; 4*u.^4-32*u.^3+88*u.^2-96*u+36];
       elseif Prob.P == 18  % Zimmermann
           x  = linspace(0,7);
           xx = [0 0.001 0.005 0.01 0.05 0.1 0.5 linspace(0,100,93)];
           cx = [x ; x ; xx ];
           cy = [2+sqrt(16-(x-3).^2) ; 2-sqrt(16-(x-3).^2) ; 14./xx];
           cy(3,1) = 100;
       elseif Prob.P == 19  % Bump
           cx = [0 u 10];
           cy = [10 0.75./u 0];
       elseif Prob.P == 22  % HGO 468:1 + constraint
           cx = 0.7+15*(u-0.6).^2;
           cy = u;
%            cx = u;
%            cy = u.^2;
       elseif Prob.P == 23
           cx = [0.2*(v-75).^2+55 ; u];
           cy = [v ; -0.2*(u-35).^2+70];
       end

   case{'glcIP_prob'}
       if Prob.P == 5  ||  Prob.P == 4
           % Polynomial used to approximate first part
           pol1  = [0.1859  -1.45   3.6023  -4.008   4.042];
           x1 = 1:0.01:2.6;
           y1 = polyval(pol1,x1);
           % Polynomial used to approximate second part
           pol2 = [0.020281 -0.489223 4.907566 -26.24658 79.13261 -128.1209 90.3686];
           x2 = 1.97:0.01:5;
           y2 = polyval(pol2,x2);
           cf1 = zeros(1,length(x1));
           for kk = 1:length(x1)
               x = [x1(kk) y1(kk)];
               cf1(kk) = eval([Prob.FUNCS.f '(x,Prob)']);
           end
           cf2 = zeros(1,length(x2));
           for kk = 1:length(x2)
               x = [x2(kk) y2(kk)];
               cf2(kk) = eval([Prob.FUNCS.f '(x,Prob)']);
           end
           if MaxFLAG
               cf1 = - cf1;
               cf2 = - cf2;
           end
           plot3(x1,y1,cf1,'k','LineWidth',2)
           plot3(x2,y2,cf2,'k','LineWidth',2)
           return
       elseif Prob.P == 6
           % Polynomial used to approximate first part
           pol1 = [-0.04675 0.73293 -4.52002 13.93376 -22.87547 19.14755 -3.23099];
           x1 = 1:0.01:4;
           y1 = polyval(pol1,x1);
           % Polynomial used to approximate second part
           pol2 = [-0.00151416 0.0578744 -0.929714237 8.1783569 -42.725748 133.0863 -229.758893 174.3955];
           x2 = 2.5175:0.01:5.8;
           y2 = polyval(pol2,x2);
           % Line to approximate last part
           k1 = (4.4533-y2(end))/(10-x2(end));
           m1 = 4.4533-k1*10;
           x3 = x2(end):0.01:10;
           y3 = k1.*x3+m1;
           cf1 = zeros(1,length(x1));
           for kk = 1:length(x1)
               x = [x1(kk) y1(kk)];
               cf1(kk) = eval([Prob.FUNCS.f '(x,Prob)']);
           end
           cf2 = zeros(1,length(x2));
           for kk = 1:length(x2)
               x = [x2(kk) y2(kk)];
               cf2(kk) = eval([Prob.FUNCS.f '(x,Prob)']);
           end
           cf3 = zeros(1,length(x3));
           for kk = 1:length(x3)
               x = [x3(kk) y3(kk)];
               cf3(kk) = eval([Prob.FUNCS.f '(x,Prob)']);
           end
           if MaxFLAG
               cf1 = - cf1;
               cf2 = - cf2;
               cf3 = - cf3;
           end
           plot3(x1,y1,cf1,'k','LineWidth',2)
           plot3(x2,y2,cf2,'k','LineWidth',2)
           plot3(x3,y3,cf3,'k','LineWidth',2)
           return
       elseif Prob.P == 10
           cx = [u ; u];
           cy = [2*u.^4-8*u.^3+8*u.^2+2; 4*u.^4-32*u.^3+88*u.^2-96*u+36];
       end
       
   case{'con_prob'}
       if Prob.P == 1  ||  Prob.P == 2  ||  Prob.P == 3
           cx1 = linspace(-10,0.99,length(u)-1);
           cx2 = linspace(1.01,10,length(u)-1);
           cx3 = linspace(-10,-1,length(u));
           cx4 = linspace(1,10,length(u));
           cx = [cx1 1 ; 1 cx2 ; cx3 ; cx4];
           cy = [(cx1-1.5)./(cx1-1) 10 ; -10 (cx2-1.5)./(cx2-1) ; -10./cx3 ; -10./cx4];
       elseif Prob.P == 5
           cx = u;
           cy = u.^2;
       elseif Prob.P == 8
           cx = -2:0.01:2;
           cy = sqrt(1-0.25*cx.^2);
           cx = [cx ; cx];
           cy = [cy ; -cy];
       elseif Prob.P == 11
           cx = 0:0.01:1;
           cy = sqrt(1-cx.^2);
       elseif Prob.P == 16
           cx = u;
           cy = 1./cx;
       end
   end  %SWITCH

   if isfield(Prob,'SCALE')
      if Prob.SCALE == 1
         cx = (cx-Prob.xL(1))/Prob.xD(1);
         cy = (cy-Prob.xL(2))/Prob.xD(2);
      end
   end
   
   for k = 1:size(cx,1)
      I = find( all([cx(k,:) ; cy(k,:) ; -cx(k,:) ; -cy(k,:)] >= ...
         repmat([Prob.x_L ; -Prob.x_U],1,size(cx,2))) );

      cxI = cx(k,I);
      cyI = cy(k,I);
      cfI = zeros(1,length(I));
      for kk = 1:length(I)
         x = [cxI(kk) cyI(kk)];
         cfI(kk) = eval([Prob.FUNCS.f '(x,Prob)']);
      end
      cfI = min(cfI,MAXW);
      
      if MaxFLAG
         cfI = - cfI;
      end
      plot3(cxI,cyI,cfI,'k','LineWidth',2)
   end
   
end
%    xlim([Prob.x_L(1) Prob.x_U(1)])
%    ylim([Prob.x_L(2) Prob.x_U(2)])
   

function plotLinearConstraints(Prob,MAXW,MaxFLAG,x_L,x_U)
% Linear constraints
A = Prob.A;

if ~isempty(A)
   [m,n]= size(A);
   b1 = Prob.b_L;
   if isempty(b1)
      b1 = -inf(m,1);
   end
   b2 = Prob.b_U;
   if isempty(b2)
      b2 = inf(m,1);
   end
   xL1 = x_L(1);
   xL2 = x_L(2);
   xU1 = x_U(1);
   xU2 = x_U(2);
else
   return
end

py = [xL2 xU2];
for k = 1:m
   if ~isinf(b2(k))
      b = b2(k);
   else
      b = b1(k);
   end
   
   a1 = A(k,1);
   a2 = A(k,2);
   
   if a1 == 0
      hLine = b/a2;
      if hLine < xU2  ||  hLine > xL2
         X = linspace(xL1,xU1);
         Y = hLine*ones(1,100);
      end
   elseif a2 == 0
      vLine = b/a1;
      if vLine < xU1  ||  vLine > xL1
         X = vLine*ones(1,100);
         Y = linspace(xL2,xU2);
      end
   else
      px = (b-a2.*py)/a1;

      if px(1) > xU1  ||  px(1) < xL1
         px(1) = min(xU1, max(xL1,px(1)) );
      end
      if px(2) > xU1  ||  px(2) < xL1
         px(2) = min(xU1, max(xL1,px(2)) );
      end
      X = linspace(px(1),px(2));
      Y = (b-a1.*X)/a2;
   end
   
   F = zeros(1,100);
   for kk = 1:100
      x = [X(kk) Y(kk)];
      F(kk) = eval([Prob.FUNCS.f '(x,Prob)']);
   end
   F = min(F,MAXW);
   if MaxFLAG
      F = -F;
   end
   plot3(X,Y,F,'k','LineWidth',2)
end

%    xlim([Prob.x_L(1) Prob.x_U(1)])
%    ylim([Prob.x_L(2) Prob.x_U(2)])


% MODIFICATION LOG:
%
% 080410 hkh Written
% 080620 nhq Added plot of linear constraints
% 080705 nhq Added glc_prob 24 constraints
% 080915 nhq Added name and number of problem to figure title
% 080916 nhq Added glc_prob 23 constraints
% 080919 nhq Added glc_prob 25 constraint

%print(27,'-depsc','testing')
%print(27,'-dpng','testing')