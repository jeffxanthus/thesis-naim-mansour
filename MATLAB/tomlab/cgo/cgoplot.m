% cgoplot.m
%
% Utility to make plots in TOMLAB /CGO
%
% function  cgoplot(task, plotData, Runs, Solver, varargin)
%
% Input:
%    task      What to do:
%       task==1   Plot contour plot with iterations
%       task==2   Mesh  plot with surfc
%       task==3   Mesh  plot with meshc
%       task==4   Plot function values for each iteration
%       task==5   Plot linear    convergence rate estimate for each iteration
%       task==6   Plot quadratic convergence rate estimate for each iteration
%       task==7   Plot weighted residual at start and after optimization
%       task==8   Plot weighted residual at start
%       task==9   Plot weighted residual after optimization
%       task==10  Plot model at start and after opt. against data
%       task==11  Plot model at start against data
%       task==12  Plot model after optimization against data
%       task==13  Plot data
%       task==14  Plot circle fitted to data point in the plane
%       task==15  Plot sampled data in global optimization
%       task==16  Plot sampled data in constrained global optimization
%       task==17  Plot contour plot and sampled data in global optimization
%       task==18  Plot contour plot and sampled data in global con optimization
%    plotData  Structure with plot information
%    Runs      Which runs to plot. Default Last entry in plotData.
%    Solver    Solver used for the run, used in global optimization
%
% Structure plotData. Used fields:
%  Prob.x_U     Upper bounds on x
%  Prob.x_L     Lower bounds on x
%  Prob.Name    Name of current test problem
%  Prob.x_opt   Set of stationary points
%  meshXcenter  Center of plot, x
%  Prob.x_0     Starting point for optimization
% task==1
%  Prob.p_dx    Matrix with search directions
%  Prob.alphaV  Matrix with steplengths
%  X_min        Lower bounds on plot domain
%  X_max        Upper bounds on plot domain
% task==2,3:
%  ixAxis       Index pointer, which search step in p_dx to use. Plotting in
%               two dimensional space spanned by x(ixAxis(1) and x(ixAxis(2))
%  viewPnt      View Point for mesh.
%  PlotPnts     Number of points in x and y direction, respectively.
% task==4:
%  F_X          Matrix with rows: [iter_no f(x)].
% task==5:
%  beta1        Linear    convergence rate. Estimated for each iter.
% task==6:
%  beta2        Quadratic convergence rate. Estimated for each iter.
% task==7,8,9:
%  Prob.LS.t  Time vector
% task==10,11,12,13:
%  Prob.LS.t  Time vector
%  Prob.LS.y  Data matrix or vector y
% task==14:
%  Prob.LS.y  Data matrix with points in the plane
%  Prob.uP      Parameters for circle problem, center point and radius
% task==15:

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1990-2006 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written May 31, 1990.    Last modified Aug 24, 2009.

function  cgoplot(task, plotData, Runs, Solver, varargin)

if nargin < 3
   Runs=length(plotData);
end

LineStyle=['r-';'b-';'g-';'c-';'m-';'r:';'b:';'g:';'c:';'m:'];
Last=Runs(1);
figs=plotData(Last).figs(:);
figs=[figs;figs(1)+1;figs(1)+2]; % For safety, if figs too short.

figure(figs(1));

Prob=plotData(Last).Prob;

if any(task == [ 1 17 18])  % Plot iterations
   Userf = Prob.FUNCS.f;  % Directly call the function, no gateway
   x_U   = Prob.x_U;
   x_L   = Prob.x_L;
   Name  = Prob.Name;
   x_opt = Prob.x_opt;
   n     = Prob.N;
   if ~isempty(Prob.optParam)
      cTol = Prob.optParam.cTol;    % Tolerance when constraint is == 0
   else
      cTol = 1E-6;
   end

   meshXcenter =  plotData(Last).meshXcenter;
   x_min =        plotData(Last).X_min;
   x_max =        plotData(Last).X_max;
   ix =           plotData(Last).ixAxis;

   if isempty(meshXcenter)
      meshXcenter=Prob.x_0;
   end

   if n==1 
      xmin = x_min(1);
      xmax = x_max(1);
      xstep = (xmax-xmin)/50;
      tc=(xmin:xstep:xmax)';
      yc=tc;
      if nargin(Userf)==1
         for i=1:length(tc)
             yc(i)=feval(Userf,tc(i));
         end
      elseif nargin(Userf)==2
         for i=1:length(tc)
             yc(i)=feval(Userf,tc(i),Prob);
         end
      else
         for i=1:length(tc)
             yc(i)=feval(Userf,tc(i),Prob, varargin{:});
         end
      end
      %title([sprintf('Plot using function %s', Userf)]); hold on;
      title(Name);
      hold on;
      plot(tc,yc);
      % Plot the optimization steps

      for k1=1:length(Runs)
         k=Runs(k1);
      
         alphak= plotData(k).alphaV;
         x0k   = plotData(k).Prob.x_0;
         dxk   = plotData(k).p_dx;
      
         if isempty(alphak)
            alphak=zeros(size(dxk,2),2);
            alphak(:,1)=1;
         end
         %% correct depending on n>=2
         %x0k = [x0k(ix(1));x0k(ix(2))];
         %dxk = [dxk(ix(1),:);dxk(ix(2),:)];
         [m,n] = size(alphak);
         x = x0k;
         for i=1:m
            j = n;
            while alphak(i,j) == 0 & j > 1
               j = j-1;
            end
            alphai = alphak(i,1:j);
            x = [x x0k*ones(1,length(alphai))+dxk(:,i)*alphai];
            x0k = x0k+dxk(:,i)*alphai(j);
         end
         yc=[];
         for i=1:size(x,2)
             yc(i)=feval(Userf,x(1,i),Prob, varargin{:});
         end
         plot(x(1,:),yc,LineStyle(k,:),x(1,:),yc,'kx');
         plot(x(1,1),yc(1),'bo');
         %plot(x(1,size(x,2)),x(2,size(x,2)),'bo');
         % Plot stationary points
         if ~isempty(x_opt)
            if size(x_opt,1)==1 % x_opt is a vector
               f_opt=feval(Userf,x_opt,Prob, varargin{:});
               x_opt=x_opt(:)'; % row
            end
            nn=1;
            if size(x_opt,1)==1
               % x_opt includes only optimal point
               plot(x_opt(1,1),f_opt,'b*');
            elseif size(x_opt,2)==nn+1
               % x_opt has nn+1 columns
               for i=1:size(x_opt,1)
                   f_opt=feval(Userf,x_opt(i,1),Prob, varargin{:});
                   if x_opt(i,nn+1)==0 % min point
                      plot(x_opt(i,ix(1)),f_opt,'b*');
                   elseif x_opt(i,nn+1)==1 % saddle point
                      plot(x_opt(i,ix(1)),f_opt,'b^');
                   elseif x_opt(i,nn+1)==2 % max point
                      plot(x_opt(i,ix(1)),f_opt,'bd');
                   end
               end
            elseif size(x_opt,2)==nn
               % x_opt has nn+1 columns
               for i=1:size(x_opt,1)
                   f_opt=feval(Userf,x_opt(i,1),Prob, varargin{:});
                   plot(x_opt(i,ix(1)),f_opt,'b*');
               end
            end
         end
      end
      hold off;
      return % Back to MAIN LOOP
   end 

   % When we reach this line, the number of variables is >=2 
   if length(ix)<2, ix = [1 2]; end
   xmin = x_min(ix(1));		
   xmax = x_max(ix(1));
   ymin = x_min(ix(2));
   ymax = x_max(ix(2));
   xstep = (xmax-xmin)/50;
   ystep = (ymax-ymin)/50;

   % Here starts "step2d.m"
   N = 40;   % Number of level curves
   logf=0;   % Do not use log transformation now
   % Init and draw contours
   xv = xmin:xstep:xmax;
   yv = ymin:ystep:ymax;
   ly = length(yv);
   lx = length(xv);
   V  = zeros(ly,lx);
   [X Y]=meshgrid(xv,yv);
   n = length(x_min);
   if n==2      % # of variables = 2
      w=[nan;nan];
   elseif n > 2 % # of variables = 2
      w=meshXcenter(:);
   end
   ix1=ix(1);
   ix2=ix(2);
   if nargin(Userf)==1
      for i=1:lx
          for j=1:ly
              w(ix1) = X(i,j);
              w(ix2) = Y(i,j);
              V(i,j) = feval(Userf, w);
          end
      end
   elseif nargin(Userf)==2
      for i=1:lx
          for j=1:ly
              w(ix1) = X(i,j);
              w(ix2) = Y(i,j);
              V(i,j) = feval(Userf, w, Prob);
          end
      end
   else
      for i=1:lx
          for j=1:ly
              w(ix1) = X(i,j);
              w(ix2) = Y(i,j);
              V(i,j) = feval(Userf, w, Prob, varargin{:});
          end
      end
   end
   if logf
      M=min(min(V));
      if M < 0
         V=V+abs(M)+1E-13;
      end
      V = log(V);
   end
   %clf;
   contour(xv,yv,V,N);
   axis([xmin ; xmax ; ymin ; ymax]);
   %title(Name);
   hold on; 

   % Plot constraints (just plot inside axis: xmax, xmin, ymax, ymin)
   w(ix(1))=xmin;
   w(ix(2))=ymin;
   [cmin,Pdim]=nlp_cF(w, Prob, varargin{:});
   me=Pdim(2);
   if ~isempty(cmin) % There are constraints  
      points=20;
      %options.TolX=sqrt(eps);
      %options.Display='off';
      xvec=linspace(xmin,xmax,points);
      yvec=linspace(ymin,ymax,points);
      Prob.FUNCS.f0='conzero';
      cMatrix=zeros(points,points,length(cmin));
      for i=1:points
          for j=1:points
             w(ix(1))=xvec(i);
             w(ix(2))=yvec(j);
             cc=nlp_cF(w, Prob, varargin{:});

             cMatrix(i,j,:)=cc; % Constraint values 
             %idx = find( abs(cc) <= cTol );
             %cMatrix(i,j,idx)=0;
             cMatrix(i,j, abs(cc) <= cTol )=0;

             %if abs(cc(index)) <= cTol
             %   cMatrix(i,j)=0;
             %else
             %   cMatrix(i,j)=cc(index); % Constraint values 
             %end
          end
      end
      for index=1:length(cmin)
         k = 0;
         if index <= me
            cCol = 'b.';
         else
            cCol = 'k.';
         end
         matrix=sign(cMatrix(:,:,index));
         %for i=1:points
         %   for ii=1:points
         %      if matrix(i,ii)==0
         %         cmat=[cmat;[index xvec(i) yvec(ii)]];
         %      elseif i>=2&ii==1
         %         if sign(matrix(i,ii))*sign(matrix(i-1,ii))==-1
         %            xzero=fzero('conzero',[xvec(i);xvec(i-1)], ...
         %                  options,yvec(ii),w,Prob,index,ix,1,...
         %                   varargin{:});
         %            cmat=[cmat;[index xzero yvec(ii)]];
         %         end
         %      elseif i==1&ii>=2
         %         if sign(matrix(i,ii))*sign(matrix(i,ii-1))==-1
         %            yzero=fzero('conzero',[yvec(ii);yvec(ii-1)], ...
         %                  options,xvec(i),w,Prob,index,ix,2,...
         %                   varargin{:});
         %            cmat=[cmat;[index xvec(i) yzero]];
         %         end
         %      elseif i>=2&ii>=2
         %         if sign(matrix(i,ii))*sign(matrix(i-1,ii))==-1
         %            xzero=fzero('conzero',[xvec(i);xvec(i-1)], ...
         %                  options,yvec(ii),w,Prob,index,ix,1,...
         %                   varargin{:});
         %            cmat=[cmat;[index xzero yvec(ii)]];
         %         end
         %         if sign(matrix(i,ii))*sign(matrix(i,ii-1))==-1
         %            yzero=fzero('conzero',[yvec(ii);yvec(ii-1)], ...
         %                  options,xvec(i),w,Prob,index,ix,2,...
         %                  varargin{:});
         %            cmat=[cmat;[index xvec(i) yzero]];
         %         end
         %      end   
         %   end
         %end
         for i=1:points
            for j=1:points
               if matrix(i,j)==0
                  k = k+1;
                  plot(xvec(i), yvec(j),cCol);
               elseif i >=2 & j==1
                  if matrix(i,j) * matrix(i-1,j) == -1
                     %xzero=fzero('conzero',[xvec(i);xvec(i-1)], ...
                     %      options,yvec(j),w,Prob,index,ix,1,...
                     %       varargin{:});

                     Prob.f0 =struct('y',yvec(j),'z',w,'index',index', ...
                                     'ix',ix,'var',1);

                     xzero=Tfzero(xvec(i-1),xvec(i), Prob,[],1E-3);

                     k = k+1;
                     plot(xzero, yvec(j),cCol);
                  end
               elseif i==1 & j>=2
                  if matrix(i,j) * matrix(i,j-1) == -1
                     %yzero=fzero('conzero',[yvec(j);yvec(j-1)], ...
                     %      options,xvec(i),w,Prob,index,ix,2,...
                     %       varargin{:});

                     Prob.f0 =struct('y',xvec(i),'z',w,'index',index', ...
                                     'ix',ix,'var',2);

                     yzero=Tfzero(yvec(j-1),yvec(j), Prob,[],1E-3);
                     k = k+1;
                     plot(xvec(i), yzero,cCol);
                  end
               elseif i>=2 & j>=2
                  if matrix(i,j) * matrix(i-1,j) == -1
                     %xzero=fzero('conzero',[xvec(i);xvec(i-1)], ...
                     %      options,yvec(j),w,Prob,index,ix,1,...
                     %       varargin{:});
                     Prob.f0 =struct('y',yvec(j),'z',w,'index',index', ...
                                     'ix',ix,'var',1);
                     xzero=Tfzero(xvec(i-1),xvec(i), Prob,[],1E-3);

                     k = k+1;
                     plot(xzero, yvec(j), cCol);
                  end
                  if matrix(i,j) * matrix(i,j-1) == -1
                     %yzero=fzero('conzero',[yvec(j);yvec(j-1)], ...
                     %      options,xvec(i),w,Prob,index,ix,2,...
                     %      varargin{:});

                     Prob.f0 =struct('y',xvec(i),'z',w,'index',index', ...
                                     'ix',ix,'var',2);

                     yzero=Tfzero(yvec(j-1),yvec(j), Prob,[],1E-3);
                     k = k+1;
                     plot(xvec(i), yzero, cCol);
                  end
               end   
            end
         end
         %fprintf('Constraint %d. Plotted %d points\n',index,k)
      end
   end
   
   % Plot upper and lower bounds
   if isempty(x_L), x_L=-Inf*ones(n,1); end
   if isempty(x_U), x_U= Inf*ones(n,1); end
   x_L=max([-1e20*ones(size(x_L(:)))';x_L(:)'])'; % Correct if -inf
   x_U=min([1e20*ones(size(x_U(:)))';x_U(:)'])';  % Correct if inf
   xbounds=[x_L(ix(1)) x_U(ix(1)) x_U(ix(1)) x_L(ix(1)) x_L(ix(1))];
   ybounds=[x_L(ix(2)) x_L(ix(2)) x_U(ix(2)) x_U(ix(2)) x_L(ix(2))];
   plot(xbounds,ybounds,'k:')

  if task == 1
   % Plot the optimization steps

   for k1=1:length(Runs)
      k=Runs(k1);
   
      alphak= plotData(k).alphaV;
      x0k   = plotData(k).Prob.x_0;
      dxk   = plotData(k).p_dx;
   
      if isempty(alphak)
         % Assume that the points in p_dx are the true x values
         if ~isempty(dxk)
            plot(dxk(ix(1),:),dxk(ix(2),:),'rs');
         end

         %alphak=zeros(size(dxk,2),2);
         %alphak(:,1)=1;
      else
         % correct depending on n>=2
         x0k = [x0k(ix(1));x0k(ix(2))];
         dxk = [dxk(ix(1),:);dxk(ix(2),:)];
         [m,n] = size(alphak);
         x = x0k;
         for i=1:m
            j = n;
            while alphak(i,j) == 0 & j > 1
               j = j-1;
            end
               alphai = alphak(i,1:j);
            x = [x x0k*ones(1,length(alphai))+dxk(:,i)*alphai];
            x0k = x0k+dxk(:,i)*alphai(j);
         end
         plot(x(1,:),x(2,:),LineStyle(k,:),x(1,:),x(2,:),'kx');
         %plot(x(1,:),x(2,:),LineStyle(k,:));
         %plot(x(1,:),x(2,:),'kx');
         plot(x(1,1),x(2,1),'bo');
         plot(x(1,size(x,2)),x(2,size(x,2)),'bo');
      end
      % Plot stationary points
      if ~isempty(x_opt)
         if size(x_opt,1)==1|size(x_opt,2)==1 % x_opt is a vector
            x_opt=x_opt(:)'; % row
         end
         nn=length(meshXcenter); % # of variables
         if size(x_opt,1)==1&size(x_opt,2)==nn
            % x_opt includes only optimal point
            plot(x_opt(1,ix(1)),x_opt(1,ix(2)),'gs');
         elseif size(x_opt,2)==nn+1
            % x_opt has nn+1 columns
            for i=1:size(x_opt,1)
                if x_opt(i,nn+1)==0 % min point
                   plot(x_opt(i,ix(1)),x_opt(i,ix(2)),'gs');
                elseif x_opt(i,nn+1)==1 % saddle point
                   plot(x_opt(i,ix(1)),x_opt(i,ix(2)),'g^');
                elseif x_opt(i,nn+1)==2 % max point
                   plot(x_opt(i,ix(1)),x_opt(i,ix(2)),'gd');
                end
            end
         elseif size(x_opt,2)==nn
            % x_opt has nn columns, assume all points are minima
            for i=1:size(x_opt,1)
                plot(x_opt(i,ix(1)),x_opt(i,ix(2)),'gs');
            end
         end
      end
   end
  end
   if task==17 | task == 18
      % Plot all sampled points in global optimization
      x_U   = Prob.x_U;
      x_L   = Prob.x_L;
      Name1 = deblank(Prob.Name);
      x_opt = Prob.x_opt;
   
      x_min = plotData(Last).X_min;
      x_max = plotData(Last).X_max;
      x_k   = plotData(Last).meshXcenter;
      ix    = plotData(Last).ixAxis;
      switch Solver
         case 'glbSolve'
            glbFile = 'glbSave.mat';
         case 'glcSolve'
            glbFile = 'glcSave.mat';
         case 'glbFast'
            glbFile = 'glbFastSave.mat';
         case {'glcFast','glcCluster'}
            glbFile = 'glcFastSave.mat';
         case {'rbfSolve','ego'}
            glbFile = 'cgoSave.mat';
         otherwise
            glbFile = '';
      end
   
      if exist(glbFile,'file')
         load(glbFile,'Name')
         if strcmpi(Name1,deblank(Name))
            if strcmpi(glbFile,'cgoSave.mat')
               load(glbFile,'O')
               plot(O(ix(1),:),O(ix(2),:),'.b');
            else
               load(glbFile,'C')
               % Transform from [0,1] to original coordinates
               C = tomsol(9, x_L, C, x_U-x_L); 
               plot(C(ix(1),:),C(ix(2),:),'.b');
            end
            s = ['Problem ' Name  ' . - sampled point'];
            hold on;
            if ~isempty(x_opt)
%if max(ix) > size(x_opt,2)
%x_opt
%ix
%   keyboard
%end
               plot(x_opt(:,ix(1)),x_opt(:,ix(2)),'or');
               s = [s ', o - known global optima'];
            end
            if ~isempty(x_k)
               plot(x_k(ix(1),:),x_k(ix(2),:),'*g');
               s = [s ', * - best point found'];
            end
            title(s);
            xlabel(['Variable ' num2str(ix(1)) ]);
            ylabel(['Variable ' num2str(ix(2)) ]);
            zoom on
            grid on
         else
            fprintf('\n')
            fprintf('Previous run was problem %s\n',Name)
            fprintf('This plot is for problem %s\n',Name1)
            fprintf('No points to plot. Optimize first\n')
         end
      else
         fprintf('CANNOT FIND the file %s. Compute global solution!\n',glbFile);
      end
   elseif ~isempty(Prob)
      title(['Problem ' Name '. Contour plot']);
      xlabel(['Variable ' num2str(ix(1)) ]);
      ylabel(['Variable ' num2str(ix(2)) ]);
   end
   hold off;

end

% --------------------------------------------------------------------

if task==2 | task==3  % Mesh  plot with surfc or meshc

   task23(task, plotData, Prob, Last, varargin{:});

% --------------------------------------------------------------------

%
% Plot function values for each iteration
%

elseif task==4

   clf;
   minF=Inf;
   for k1=1:length(Runs)
       k=Runs(k1);
       minF=min(minF,min(plotData(k).F_X(:,2)));
   end
   if minF > 0, minF=0; else minF=abs(minF)+1E-15; end
   if minF == 0
      title([ Prob.Name '. f(x) each iteration']);
   else
      title([ Prob.Name ...
              '. Added ' num2str(minF) ' + 1E-15 to f(x)']);
   end
   ylabel('Function values: log10(f(x))');
   xlabel('Iteration number');
   hold on;
   for k1=1:length(Runs)
       k=Runs(k1);
       plot(plotData(k).F_X(:,1),log10(minF+plotData(k).F_X(:,2)),...
           [LineStyle(k,1) 'x']); 
   end
   hold off;

% --------------------------------------------------------------------

%
% Plot convergence rate, both linear and quadratic
%

elseif task==5
   clf;
   title([ Prob.Name '. Linear convergence rate']);
   ylabel('log10(Convergence rate)');
   xlabel('Iteration number');
   hold on;
   for k1=1:length(Runs)
       k=Runs(k1);
       plot(1:length(plotData(k).beta1),log10(plotData(k).beta1),...
           [LineStyle(k,1) 'x']); 
   end
   hold off;

elseif task==6
   clf;
   title([ Prob.Name '. Quadratic convergence rate']);
   ylabel('log10(Convergence rate)');
   xlabel('Iteration number');
   hold on;
   for k1=1:length(Runs)
       k=Runs(k1);
       plot(1:length(plotData(k).beta2),log10(plotData(k).beta2),...
           [LineStyle(k,1) 'x']); 
   end
   hold off;

% --------------------------------------------------------------------

% Plot residual for least squares problems, 
% The residual both at start and after optimization 
%
elseif task==7 | task==8 | task==9

   clf;
   Name=Prob.Name;
   if task==7
      title(['Residual ' Name ...
             '. Color: ' LineStyle(1:length(Runs),1)' ...
             '. At start "o". After optimization "*"']);
   elseif task==8
       title(['Residual at start. Problem ' Name ... 
            '. Color: ' LineStyle(1:length(Runs),1)' ...
            ]);
   elseif task==9
       title(['Residual after optimization. Problem ' Name ...
            '. Color: ' LineStyle(1:length(Runs),1)' ...
            ]);
   end

   hold on;
   t=Prob.LS.t(:);
   if isempty(t)
      xlabel('Point number');
   else
      xlabel('Time t');
   end

   for k1=1:length(Runs)
       k=Runs(k1);
       x_0=plotData(k).Prob.x_0;
       x_k=plotData(k).meshXcenter;

       if task==7 | task==9
          r_0=nlp_r(x_0, Prob, varargin{:});
       else
          r_0=nlp_r(x_k, Prob, varargin{:});
       end
      

       if task==7
          if any(x_0~=x_k)
             r_k = nlp_r(x_k, Prob, varargin{:});
             if isempty(t)
                plot(r_k,[LineStyle(k,1) '*']);
             else
                plot(t,r_k,[LineStyle(k,1) '*']);
             end
          end
       end
       if isempty(t)
          xlabel('Point number');
          if task==9
             plot(r_0,[LineStyle(k,1) '*']);
          else
             plot(r_0,[LineStyle(k,1) 'o']);
          end
       else
          xlabel('Time t');
          if task==9
             plot(t,r_0,[LineStyle(k,1) '*']);
          else
             plot(t,r_0,[LineStyle(k,1) 'o']);
          end
       end
   end
end

% --------------------------------------------------------------------

%
% Plot model and data for least squares (LS) problems --------------
%
% Plot how the LS model fits to data using:
% 1. the starting values. 
% 2. the optimal values after optimization
%
% The routine yModel is used to compute the model values
%
% The model is plotted with higher resolution. Default N=201 points
%

if task==10 | task==11 | task==12 | task==13 
   Name=Prob.Name;
   P=Prob.P;

   x_0=Prob.x_0;
   x_k=plotData(Last).meshXcenter;

   y=Prob.LS.y;
   t=Prob.LS.t(:);
   tM=[];

   m=size(y,1);

   if task==10 | task==12
      if isempty(t)
         [yMod tM]=yModel('nlp_r',x_k,Prob,[1, length(y)]);
         if length(tM)~=length(yMod)
            if length(yMod)==length(y)
               tM=1:length(y);
            else
               tM=[]; 
            end
         end
      else
         [yMod tM]=yModel('nlp_r',x_k,Prob,[ t(1) t(m)]);
      end
   end

   if task==10 | task==11
      if any(x_k~=x_0) | task==11
         if isempty(t)
            if isempty(tM)
               [yM1 tM]=yModel('nlp_r',x_0,Prob,[1, length(y)]);
            else
               yM1=yModel('nlp_r',x_0,Prob,tM);
            end
         elseif isempty(tM)
            [yM1 tM]=yModel('nlp_r',x_0,Prob,[t(1) t(m)]);
         else
            yM1=yModel('nlp_r',x_0,Prob,tM);
         end
      else
         yM1=yMod;
         yMod=[];
      end
   end
   if task==11
      yMod=[];
   end
   if task==12
      yM1=yMod;
      yMod=[];
   end
   if task==13
      yMod=[];
      yM1=[];
      tM=[];
   end

   N=length(y);
   m=round(N/2);
   figure(figs(1));
   plotymod(['Problem # ' num2str(P) '. ' Name], t, y, tM, yMod, yM1);
   if m > 6 & task ~= 13
      figure(figs(2));
      plotymod(['Problem # ' num2str(P) '. ' Name '; 1st part'], ...
               t, y, tM, yMod, yM1, 1, m);
   
      figure(figs(3));
      plotymod(['Problem # ' num2str(P) '. ' Name '; 2nd part'], ...
               t, y, tM, yMod, yM1, m+1, N);
   end
   hold off;
   
end

% --------------------------------------------------------------------

% Plot approximating circle for the the circle problem
%

if task==14
   uP=Prob.uP;
   y=Prob.LS.y;
   m=size(y,1);
   if size(y,2) == 1
      m=m/2;
      y=reshape(y,m,2);
   end

   CentCi=uP(1:2);  % Center for theoretical circle
   radius=uP(3);    % Radius for theoretical circle
   clf;
   title(['Approx circle. Color: ' LineStyle(1:length(Runs),1)' ...
          ' (' num2str(m) ' points). Theoretical (red/dashed) ']); hold on;
   xlabel('x'); ylabel('y');
   plot(CentCi(1),CentCi(2),'or'); 
   tt=[0:0.05:2]'*pi+0.1;
   xo=cos(tt)*radius(1)+CentCi(1);
   yo=sin(tt)*radius(1)+CentCi(2);
   plot(xo,yo,'--r');
   plot(y(:,1),y(:,2),'+');

   for k1=1:length(Runs)
       k=Runs(k1);
       x_k=plotData(k).meshXcenter;
       xc=cos(tt)*x_k(1)+x_k(2);
       yc=sin(tt)*x_k(1)+x_k(3);
       plot(xc,yc,LineStyle(k,:)); 
   end
   
   hold off;
end
if task==15 | task==16

   % Plot all sampled points in global optimization
   x_U   = Prob.x_U;
   x_L   = Prob.x_L;
   Name1 = deblank(Prob.Name);
   x_opt = Prob.x_opt;

   x_min = plotData(Last).X_min;
   x_max = plotData(Last).X_max;
   ix    = plotData(Last).ixAxis;
   x_k   = plotData(Last).meshXcenter;

   switch Solver
      case 'glbSolve'
         glbFile = 'glbSave.mat';
      case 'glcSolve'
         glbFile = 'glcSave.mat';
      case 'glbFast'
         glbFile = 'glbFastSave.mat';
      case {'glcFast','glcCluster'}
         glbFile = 'glcFastSave.mat';
      case {'rbfSolve','ego'}
         glbFile = 'cgoSave.mat';
      otherwise
         glbFile = '';
   end

   if exist(glbFile,'file')
      load(glbFile,'Name')
      if strcmpi(Name1,deblank(Name))
         if strcmpi(glbFile,'cgoSave.mat')
            load(glbFile,'O')
            plot(O(ix(1),:),O(ix(2),:),'.b');
         else
            load(glbFile,'C')
            % Transform from [0,1] to original coordinates
            C = tomsol(9, x_L, C, x_U-x_L); 
            plot(C(ix(1),:),C(ix(2),:),'.b');
         end

         %clf;
         hold on;
         axis([x_min(ix(1)); x_max(ix(1)); x_min(ix(2)); x_max(ix(2));]);
         s = ['Problem ' Name  ' . - sampled point'];
         if ~isempty(x_opt)
            plot(x_opt(:,ix(1)),x_opt(:,ix(2)),'or');
            s = [s ', o - known global optima'];
         end
         if ~isempty(x_k)
            plot(x_k(ix(1),:),x_k(ix(2),:),'*g');
            s = [s ', * - best point found'];
         end
         title(s);
         xlabel(['Variable ' num2str(ix(1)) ]);
         ylabel(['Variable ' num2str(ix(2)) ]);
         zoom on
         grid on
      else
         fprintf('\n')
         fprintf('Previous run was problem %s\n',Name)
         fprintf('This plot is for problem %s\n',Name1)
         fprintf('No points to plot. Optimize first\n')
      end
   else
      fprintf('CANNOT FIND the file %s. Compute global solution!\n',glbFile);
   end

end

% --------------------------------------------------------------------
% SUB FUNCTIONS: plotymod, yModel

%
% function plotymod(Title, t, y, tM, yMod, yMod0,i,k);
%
% Plot nonlinear least squares model f(x,t) and data y(t). t may be empty.
%
% INPUT:
% Title    Plot title
% t        Time t or empty
% y        Function values y(t) at t or N-vector y
% tM       Time t for yMod0 and yMod vector
% YMod     Model values f(t,x*) at t,x* (optimal point of optimization) or 
%          Model N-vector f(x*) at optimal point x*
% YMod0    Model values f(t,x0) at t,x0 (start of optimization) or 
%          Model N-vector f(x0) at start x0
% i,k      Plot from point i to point k. Default all points
%
% OUTPUT:
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1990-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.6.0$
% Written Sept 18, 1995.   Last modified Jan 17, 2005.

function plotymod(Title, t, y, tM, yMod, yMod0,i,k)

if nargin < 8
   k=[];
   if nargin < 7
      i=1;
      if nargin < 6
         yMod0=[];
         if nargin < 5
            yMod=[];   
            if nargin < 4
               tM=[];   
            end
         end
      end
   end
end

t=t(:);
N=length(y);
if isempty(k), k=N; end

clf;
title(Title); hold on;
if isempty(t)
   ylab='y';
   xlabel('Point number');
   plot(i:k,y(i:k),'*');
   if ~isempty(tM)
      i1=min(find(i <= tM));
      k1=max(find(k >= tM));
   end
   if ~isempty(yMod)
      ylab=[ylab ', f(x*) [solid]'];
      plot(tM(i1:k1),yMod(i1:k1));
   end
   if ~isempty(yMod0)
      ylab=[ylab ', f(x0) [dashdot]'];
      plot(tM(i1:k1),yMod0(i1:k1),'-.');
   end
   ylabel(ylab);
else
   ylab='y(t)';
   xlabel('Time t');
   plot(t(i:k),y(i:k),'*');
   if ~isempty(tM)
      i1=min(find(t(i) <= tM));
      k1=max(find(t(k) >= tM));
   end
   if ~isempty(yMod)
      ylab=[ylab ', f(x*,t) [solid]'];
      plot(tM(i1:k1),yMod(i1:k1));
   end
   if ~isempty(yMod0)
      ylab=[ylab ', f(x0,t) [dashdot]'];
      plot(tM(i1:k1),yMod0(i1:k1),'-.');
   end
   ylabel(ylab);
end

% --------------------------------------------------------------------

%		yModel.m
%
% Evaluate model function for grid of points
%
% function [y, t] = yModel (Fname, x, Prob, tM);
%
% INPUT:  
%         Fname  Name of function which evaluates the residuals, e.g. ls_r.m
%         x      Parameter vector x
%         Prob   Problem structure
%         tM     Time values t to evaluate model f(t,x) for.
%                Default N=200 points are used (+ 1 end point)
%                if isempty(tM)>,  use interval [t(1),t(m)] with N+1 points
%                if length(tM)==1, use interval [0,tM(1)] with N+1 points
%                if length(tM)==2, use interval [tM(1),tM(2)] with N+1 points
%                if length(tM)==3, use interval [tM(1):tM(2):tM(3)]
%                if length(tM)>=4, use given points tM.
% OUTPUT: 
%         y      y=f(t,x)
%         t      Time vector t
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1990-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.6.0$
% Written Feb 11, 1998.    Last modified Jan 17, 2005.

function [y, t] = yModel (Fname, x, Prob, tM, varargin)

N=200;
if nargin < 4
   tM=[];
end

tM=tM(:);

if length(tM)==1     % Assume intervall [0, abs(tM)]
   t=(0:abs(tM)/N:tM)';
elseif length(tM)==2 % Assume interval [tM(1), tM(2)]
   t=(min(tM):(max(tM)-min(tM))/N:max(tM));
elseif length(tM)==3 % Assume interval [tM(1):tM(2):tM(3)] == [min:step:max]
   t=(tM(1):tM(2):tM(3));
elseif length(tM) > 4
   t=tM;
end

yx=Prob.LS.y;
tx=Prob.LS.t(:);
if isempty(tx)
   m=length(yx);
   if m==0
      % Can not do any plot when data vector is empty
      return
   end
   tx=(1:m)';
end

% uP=Prob.uP; % NOT USED

global probType
global LS_x
LS_x=[];

if probType ~=7
   if isempty(tM)
      %t=[ tx(1):(tx(length(tx))-tx(1))/N:tx(length(tx))];
      t=linspace(tx(1),tx(length(tx)),N);
   end

   %y=zeros(length(t),1);
   
   %Prob.LS.yUse=1;

   %Prob.LS.y=y;
   %Prob.LS.t=t(:);


   if strcmp(Fname,'nlp_r');
      r=nlp_r(x, Prob, varargin{:});
   else
      r=feval(Fname, x, Prob, varargin{:});
   end
   % safeguard against short series
   %if length(y) ~= length(t) & ~isempty(tM), t=tM(:); end  

   y=r+yx;
   y=spline(tx,y,t);
else
   if strcmp(Fname,'nlp_r');
      r=nlp_r(x, Prob, varargin{:});
   else
      r=feval(Fname, x, Prob, varargin{:});
   end
   y=r+yx;
   t=tx;
end


%if 0
%   plot(t,y);
%   pause
%end


% ====================================================================
function task23(task, plotData, Prob, Last, varargin)
% ====================================================================

if task==2 | task==3  % Mesh  plot with surfc or meshc

   omega    = [plotData(Last).X_max(:)';plotData(Last).X_min(:)'];
   Name     = Prob.Name;
   viewPnt  = plotData(Last).viewPnt;
   ix       = plotData(Last).ixAxis;
   nof_lines= plotData(Last).PlotPnts;
   w        = plotData(Last).meshXcenter;
   Userf    = Prob.FUNCS.f;

   if length(viewPnt)   < 2, viewPnt = [-37.5 30]; end
   if length(nof_lines) < 2, nof_lines = [60 60]; end
   if length(ix)        < 2, ix = [1 2]; end
   if length(w)         < 2, w  = [1 1]; end

   w=w(:);

   x = omega(1,1):((omega(2,1) - omega(1,1)) / (nof_lines(1) - 1)):omega(2,1);
   y = omega(1,2):((omega(2,2) - omega(1,2)) / (nof_lines(2) - 1)):omega(2,2);
   z = zeros(length(x),length(y));

   ix1=ix(1);
   ix2=ix(2);

   if nargin(Userf)==1
      for i=1:nof_lines(1)
         for j=1:nof_lines(2)
            w(ix1) = x(i);
            w(ix2) = y(j);
            z(i,j) = feval(Userf, w);
         end
      end
   elseif nargin(Userf)==2
      for i=1:nof_lines(1)
         for j=1:nof_lines(2)
            w(ix1) = x(i);
            w(ix2) = y(j);
            z(i,j) = feval(Userf, w, Prob);
         end
      end
   else
      for i=1:nof_lines(1)
         for j=1:nof_lines(2)
            w(ix1) = x(i);
            w(ix2) = y(j);
            z(i,j) = feval(Userf, w, Prob, varargin{:});
         end
      end
   end

   if task==2
      if isempty(Name)
         title(['Surfc plot - View point ' num2str(viewPnt(1)) ...
                ' and ' num2str(viewPnt(2))]);
      else
         title(['Problem ' Name '. Surfc plot - View point ' ...
                num2str(viewPnt(1)) ' and ' num2str(viewPnt(2))]);
      end
      view(viewPnt(1),viewPnt(2))
      surfc(x,y,z');
      %surf(x,y,z');
      view(viewPnt(1),viewPnt(2))
      if isempty(Name)
         title(['Surfc plot - View point ' num2str(viewPnt(1)) ...
                ' and ' num2str(viewPnt(2))]);
      else
         title(['Problem ' Name '. Surfc plot - View point ' ...
                num2str(viewPnt(1)) ' and ' num2str(viewPnt(2))]);
      end
   elseif task==3
      meshc(x,y,z');
      %mesh(x,y,z');
      view(viewPnt(1),viewPnt(2))
      %colormenu;
      if isempty(Name)
         title(['Meshc plot - View point ' num2str(viewPnt(1)) ...
                ' and ' num2str(viewPnt(2))]);
      else
         title(['Problem ' Name '. Meshc plot - View point ' ...
                num2str(viewPnt(1)) ' and ' num2str(viewPnt(2))]);
      end
   end
   xlabel(['Variable ' num2str(ix(1)) ]);
   ylabel(['Variable ' num2str(ix(2)) ]);
end

% MODIFICATION LOG:
%
% 981011  hkh  Setting default on Runs if argument is missing
%              Skipping Prob as input, instead Prob is in plotData struct
% 981018  hkh  Get y from Prob.LS.y, not Prob.y. Same for t.
% 981022  hkh  New option to plot only the data, without model, Put as #13
%              Safeguard against short series in yModel.
%              Improve handling of plot of initial model and data
% 981118  hkh  Use new flag Prob.LS.yUse
% 981127  hkh  Doing (:) on a number of places where time t is set.
% 981128  hkh  Setting LS_x empty, other wrong residual is returned, because
%              we manipulate the t vector.
% 981129  hkh  Title Quadratic (o) changed to only Quadratic
% 981208  hkh  Use spline to get model curve instead of calling with dense
%              grid of points with y=0, as the later may easily cause errors.
% 011030  hkh  Generalize plotting of global optimization sampled points to
%              all mat-files
% 011111  hkh  Add rbfSolve and glcCluster files
% 020423  hkh  cgoSave.mat now used for rbfSolve and ego
% 020702  hkh  Use strcmpi, not strcmp, for Name check
% 050117  med  mlint revision
% 060814  med  FUNCS used for callbacks instead
% 090821  hkh  mlint minor cleanup
