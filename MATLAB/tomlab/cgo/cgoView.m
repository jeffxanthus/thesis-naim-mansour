% function OutR = cgoView(R)
%
% Print a summary table of the results from running the CGO solvers
% using the RBF interpolation algorithms
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2007-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Sept 22, 2007. Last modified August 24, 2009.

function OutR = cgoView(R)

if nargin < 1
   error('cgoView needs one parameter, the Result structure');
end

OutR = [];
x_L   = R.Prob.x_L;
x_U   = R.Prob.x_U;
fTol  = R.Prob.optParam.eps_f;     % Relative convergence tolerance in f(x)
epsX  = max(1E-3,fTol);
x_opt = R.Prob.x_opt;

Its   = R.CGO.Its;
Info  = R.CGO.WarmStartInfo;
fGoal = R.CGO.fGoal;

X       = Info.X;
F       = Info.F;
% F_m     = Info.F_m; % - Not used

fMin    = Its.fMin;
fNew    = Its.fNew;
xInf    = Its.xInf;
fInf    = Its.fInf;
fnStar  = Its.fnStar;
nFunc   = Its.nFunc;
NN      = Its.n;
minSn   = Its.minSn;
minSny  = Its.minSn_y;
snNew   = Its.snNew;
modN    = Its.modN;
RelErr  = Its.RelErr;
% doM     = Its.distxNew2xMin;  % Not used for printing right now
% SoO     = Its.distxNew2minSn; % Not used for printing right now
maxN    = max(modN);

% Find max Lipschitz constant
Dist    = tomsol(30,X,X);
LipMx   = 0;
for i = 1:size(X,2)-1
    LipMx = max(LipMx,max(abs(F(i)-F(i+1:end))./Dist(i+1:end,i)));
end


N = length(Its.fNew);

fprintf('fGoal %12.6f fBest %12.6f ',fGoal, fMin(end));
if fGoal ~= 0
   fprintf('RelErr %7.3f%%',100*(fMin(end)-fGoal)/fGoal);
else
   fprintf('RelErr %7.3f%%',100*(fMin(end)-fGoal));
end
fprintf('\n')
fprintf('\n')
fprintf('Fun  N    fMin     fNew oB  -fnStar   -minSn     snNew Ac/Pred ');
if ~isempty(x_opt)
   fprintf('   dXO    snOptf ');
   %fprintf('|snOptg|  ');
   dXO    = Its.dXO;
   snOptf = Its.snOptf;
   %snOptg = Its.snOptg; % Not used for printing right now
   %snOptE = Its.snOptE; % Not used for printing
end
fprintf('   fInf oB  doInf ');
fprintf('  dInf-1   dInf--  dInf-X   LipMax    LipMx  LipNear ');
%HKH
%fprintf('     doM      SoO');
fprintf('\n')
J = nFunc(1)-1;
for i=1:N
    k = nFunc(i);
    if i == 1
       J = J+1;
       Update = 1;
    else
       if NN(i) > NN(i-1)
          J = J+1;
          Update = 1;
       else
          Update = 0;
       end
    end
    %fprintf('%3d ',Its.Iter(i))
    if Update
       fprintf('%3d ',k);
    else
       fprintf('%3dX',k);
    end
    fprintf('%2d ',modN(i))
    fprintf('%7.2f ',fMin(i))
    fprintf('%8.2g ',fNew(i))
    onB = nOnBound(X(:,J),x_L,x_U,epsX);
    fprintf('%2d ',onB)
    fprintf('%8.2g ',-fnStar(i))
    fprintf('%8.2g ',-minSn(i))
    fprintf('%9.2g ',snNew(i))
    if modN(i) == maxN
       AcPred = (snNew(i)-fNew(i))/snNew(i);
    else
       AcPred = (snNew(i)-fNew(i))/(snNew(i)-fnStar(i));
    end
    if AcPred > 0
        fprintf('%6.2f%% ',100*AcPred);
    else
        fprintf('%7.2f ',AcPred);
    end
    if ~isempty(x_opt)
       fprintf('%6.2f ',dXO(i));
       fprintf('%9.2f ',snOptf(i));
       %fprintf('%8.2g ',norm(snOptg(:,i)));
    end
    if fInf(i) < -1E6
       fprintf('%7.2f ',-inf);
    else
       fprintf('%7.2f ',fInf(i))
    end
    onBInf = nOnBound(xInf(:,i),x_L,x_U,epsX);
    fprintf('%2d ',onBInf)
    fprintf('%6.2f ',norm(xInf(:,i)-minSny(:,i)));
    if i > 1
       fprintf('%8.2g ',norm(xInf(:,i)-xInf(:,i-1)));
    else
       fprintf('         ');
    end
    if i > 1
       doI = tomsol(30,xInf(:,i),xInf(:,1:i-1));
       doInf = min(doI);
       %doInf = norm(xInf(:,i)-xInf(:,1));
       %for j=2:i-1
       %    doI = norm(xInf(:,i)-xInf(:,j));
       %    if doI < doInf
       %       doInf = doI;
       %    end
       %end
       fprintf('%8.2g ',doInf);
    %elseif i > 1
    %   fprintf('%8.2g ',norm(xInf(:,i)-xInf(:,i-1)));
    else
       fprintf('         ');
    end
    doI = tomsol(30,xInf(:,i),X(:,1:J-1));
    doXInf = min(doI);
    %doXInf = norm(xInf(:,i)-X(:,1));
    %for j=2:k-1
    %    doI = norm(xInf(:,i)-X(:,j));
    %    if doI < doXInf
    %       doXInf = doI;
    %    if doXInf == 0, break; end
    %    end
    %end
    fprintf('%7.2f ',doXInf);

    Lip = abs(F(J)-F(1:J-1))./tomsol(30,X(:,J),X(:,1:J-1))';
    Lmax = max(Lip);
    LipMx = max(LipMx,Lmax);
    %Lmax = abs(F(j)-F(1))/norm(X(:,j)-X(:,1));
    %for j=2:k-1
    %    Lip = abs(F(k)-F(j))/norm(X(:,k)-X(:,j));
    %    if Lip > Lmax
    %       Lmax = Lip;
    %    end
    %end
    fprintf('%8.2g ',Lmax);
    fprintf('%8.2g ',LipMx);

    D    = tomsol(30,X(:,J),X(:,1:J-1));
    [Dmin,ix] = min(D); 
    %Dmin = norm(X(:,k)-X(:,1));
    %ix = 1;
    %for j=2:k-1
    %    D = norm(X(:,k)-X(:,j));
    %    if D < Dmin
    %       Dmin = D;
    %       ix   = j;
    %    end
    %end
    %L = abs(F(k)-F(ix))/norm(X(:,k)-X(:,ix));
    L = abs(F(J)-F(ix))/Dmin;
    fprintf('%8.2g ',L);
    %fprintf('%8.2g ',doM(i));
    %fprintf('%8.2g ',SoO(i));

    fprintf('\n')
end
return

% Not used code
%if fGoal < 0
%   f0 = -10^(ceil(log10(abs(fGoal))));
%else
%   f0 = 10^(floor(log10(fGoal)));
%end
%v1 = log10(fMin-f0);
%v2 = log10(fNew-f0);
%v0 = log10(fGoal-f0);
%v3 = log10(minSn-f0);
%v4 = log10(snNew-f0);
%figure(1)
%plot(nFunc,v1,'x')
%hold on;
%title (['log10 of fMin(*), fNew(o), fGoal - offset ' num2str(f0) ])
%plot(nFunc,ones(1,N)*v0)
%plot(nFunc,v2,'o')
%plot(nFunc,v3,'b')
%%plot(nFunc,v4,'+')
%%%hold on;
%%plot(fNew,'o')
%figure(2)
%plot(nFunc,log10(RelErr))
%title('log10 of Relative Error')


% MODIFICATION LOG:
%
% 080410 hkh Written 
% 090824 hkh Minor mlint cleanup
