%function X = getMaxiMinLHDnorm2(dim,nSample)
%
% Optimal 2-Norm Maximin LHDs. All designs are taken
% from  http://www.spacefillingdesigns.nl/
%
%   Available Designs     dim     nSample
% ----------------------------------------
%   Maximin L2-norm       2-4     2-300
%                        5-10     2-100

function X = getMaxiMinLHDnorm2(dim,nSample)

if dim > 4
   load MaxiMinLHD5_10.mat Xs
   n = 4;
else
   load MaxiMinLHD2_4.mat Xs
   n = 1;
end
X = Xs(dim-n,nSample-1).X;
clear Xs

% MODIFICATION LOG:
%
% 080707  nhq  Written
% 080717  hkh  Changed to mat file format

