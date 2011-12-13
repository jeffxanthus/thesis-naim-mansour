% function X = getLHDnormAE(dim,nSample)
%
% Optimal Audze-Eglais LHDs. All designs are taken
% from  http://www.spacefillingdesigns.nl/
%
%   Available Designs     dim     nSample
% ----------------------------------------
%   Audze-Eglais (AE)    2-10     2-100

function X = getLHDnormAE(dim,nSample)

load LHDnormAE.mat

X = Xs(dim,nSample).X;

clear Xs


% MODIFICATION LOG:
%
% 080707  nhq  Written
% 080717  hkh  Changed to mat file format
