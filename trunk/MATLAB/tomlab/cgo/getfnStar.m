% --------------------------------------------------------------
function [fnStar, ucOK, nmax] = getfnStar(modN, fnStar, n, nmax, nFunc, nInit, N, ...
          DeltaRule, REPLACE, F, F_m, min_sn, fMin, ...
          f_Low, eps_sn, fStarRule)
% --------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Oct 6, 2009. Last modified Oct 6, 2009.

% NOTE! Hard coded constants determining current cycle strategy.
% Should be further tested and developed to see which one is best.
% The old rbfSolve

OldStrat = 0; % If false, use a new Idea 1 median strategy
% Possibly amplify fDiff to generate wider range of target values
AmpFac   = 1.0;
ucOK     = 1;

if modN == 0
   nmax = n;
else
   if DeltaRule
      % Current idea used:
      nmax = min(n,max(2,nmax-floor((nFunc-nInit)/N)));
   else
      % Always use ALL points
      nmax = n;
   end
end

if OldStrat
   F_sort = sort(F_m);
   max_F  = F_sort(nmax);
else
   if REPLACE > 1 & nmax > n/2
      max_F = median(F);
   elseif REPLACE == 1 & nmax > n/2
      max_F = median(F);
   elseif REPLACE == 0
      F_sort = sort(F);
      max_F = F_sort(nmax);
   else
      F_sort = sort(F_m);
      max_F  = F_sort(nmax);
   end
end

fOld = fnStar;

%if f_Low > -1E299
%   fDiff = max(0,min(max_F-min_sn,min_sn-f_Low));
%else
%   fDiff = max_F-min_sn;
%end
%fDiff = max(fDiff,max_F-fMin);
% Use fMin instead of min_sn, because surface way be wild
if f_Low > -1E299
   fDiff = max(0,min(max_F-fMin,min_sn-f_Low));
else
   fDiff = max_F-fMin;
end
%HKH NOTE **********************
% Special safe guard if not FEASIBLE, max(F) < min_sn
if fDiff <= 0
   %HKH New
   fDiff = max(0,max(F)-fMin);
   if fDiff <= 0
      fDiff = max(1E-4,max(F)-min(F));
   end
end

fDiff = AmpFac * fDiff;

if fStarRule == 1
   fnStar = min_sn - ( (N-modN)/N )^2*fDiff;
   % Check that the above is exactly the following:
   %fnStar = min_sn - ( (mod(N-(nFunc-nInit),N+1))/N )^2*(max_F-min_sn);
elseif fStarRule == 2
   fnStar = min_sn - ( (N-modN)/N )*fDiff;
elseif fStarRule == 3
   if modN == 0
      fnStar = min_sn;
   elseif abs(min_sn) < 1E-12
      fnStar = min_sn - ( (N-modN)/N );
   else
      fnStar = min_sn - ( (N-modN)/N )*abs(min_sn);
   end
else
   %if f_Low > -1E299
   %   fDiff = max(0,min(max_F-min_sn,min_sn-f_Low));
   %else
   %   fDiff = max_F-min_sn;
   %end
   % Use fMin instead of min_sn, because surface way be wild
   if f_Low > -1E299
      fDiff = max(0,min(max_F-fMin,min_sn-f_Low));
   else
      fDiff = max_F-fMin;
   end
   if fDiff <= 0
      %HKH New
      fDiff = max(0,max(F)-fMin);
      if fDiff <= 0
         fDiff = max(1E-4,max(F)-min(F));
      end
   end
   switch modN
      case 0
         fnStar = min_sn - fDiff;
      case 1
         fnStar = min_sn - 0.5*fDiff;
      case 2
         fnStar = min_sn - 1E-1*fDiff;
      case 3
         fnStar = min_sn - 1E-2*fDiff;
      case 4
         fnStar = min_sn - 1E-4*fDiff;
      case 5
         fnStar = min_sn;
   end
end
% fprintf('Median=max_F %20.10f \n',max_F)

if REPLACE == 1 & ...
      abs(fOld - fnStar) < 1E-4*max(abs(fOld),abs(fnStar)) & (fStarRule < 3)
   disp('ALARM. fnStar nearly equal to fnStar in last step');
   fprintf('fnStar %18.12f ',fnStar);
   fprintf('fMin %18.12f ',fMin);
   fprintf('fDiff %18.12f ',fDiff);
   fprintf('max_F %18.12f ',max_F);
   fprintf('min_sn %18.12f ',min_sn);
   fprintf('\n');
   nmax = nmax - 1;
   % Create fnStar without REPLACE
   F_sort = sort(F);
   max_F = F_sort(nmax);
   fnStar = min_sn - ( (N-modN)/N )^2*(max_F-min_sn);
   fprintf('New fnStar %18.12f ',fnStar);
   fprintf('\n');
   % Check that the above is exactly the following:
   % fnStar = min_sn - ( (mod(N-(nFunc-nInit),N+1))/N )^2*(max_F-min_sn);
end

if modN == N
   % Unconstrained cycle step. But min_sn must be sufficiently lower

   if fMin == 0
      maxF = max(F);
      if min_sn >= -eps_sn*min(1,maxF)
         if maxF == 0
            fnStar = -0.01;
         else
            fnStar = -0.01*min(1,maxF);
         end
         ucOK = 0;
      end
   elseif min_sn >= fMin-eps_sn*abs(fMin)
      fnStar = min_sn - 1E-2*abs(fMin);
      ucOK = 0;
   end
end

% MODIFICATION LOG
%
% 091002  hkh  First version written from code in rbfSolve
