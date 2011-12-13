% lgo1_prob: Defines global optimization problems. One dimension.
%
% The problems are defined by:
% Prob = probInit('lgo1_prob', i); i = 1,...,27
%
% function [probList, Prob] = lgo1_prob(P, ask, Prob);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Reference:
% Pintér, J.D., Bagirov, A., and Zhang, J. (2003) An Illustrated Collection of
% Global Optimization Test Problems. Research Report, Pintér Consulting Services,
% Inc. Halifax, NS, Canada; and CIAO-ITMS, University of Ballarat, Ballarat, Vic.,
% Australia.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 2004.   Last modified Jun 3, 2008.

function [probList, Prob] = lgo1_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'P1'...
    ,'P2'...
    ,'P3'...
    ,'P4'...
    ,'P5'...
    ,'T1'...
    ,'T2'...
    ,'T3'...
    ,'T4'...
    ,'T5'...
    ,'T6'...
    ,'T7'...
    ,'T8'...
    ,'T9'...
    ,'T10'...
    ,'T11'...
    ,'T12'...
    ,'T13'...
    ,'T14'...
    ,'T15'...
    ,'T16'...
    ,'T17'...
    ,'T18'...
    ,'T19'...
    ,'T20'...
    ,'T21'...
    ,'T22'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

if P == 1
    Name = 'P1';
    x_L = -5;
    x_U = 7;
    x_0 = 0;
    x_opt = 3.69752;
    f_opt = -145.063;
    f_Low=-200;
    x_min = -5;
    x_max = 7;
    nGlobal = 1;
    nLocal  = 1;
elseif P == 2
    Name = 'P2';
    x_L = -3;
    x_U = 8;
    x_0 = 0;
    x_opt = 6.32565;
    f_opt = -443.673;
    f_Low=-500;
    x_min = -3;
    x_max = 8;
    nGlobal = 1;
    nLocal  = 2;
elseif P == 3
    Name = 'P3';
    x_L = -2;
    x_U = 2.5;
    x_0 = 0;
    x_opt = 1.75767;
    f_opt = -0.686072;
    f_Low=-500;
    x_min = -2;
    x_max = 2.5;
    nGlobal = 1;
    nLocal  = 2;
elseif P == 4
    Name = 'P4';
    x_L = 0;
    x_U = 3;
    x_0 = 0;
    x_opt = 2;
    f_opt = -89;
    f_Low=-500;
    x_min = 0;
    x_max = 3;
    nGlobal = 1;
    nLocal  = 2;
elseif P == 5
    Name = 'P5';
    x_L = 0;
    x_U = 10;
    x_0 = 0;
    x_opt = [1;3;5;7];
    f_opt = 0;
    f_Low=-500;
    x_min = 0;
    x_max = 10;
    nGlobal = 4;
    nLocal  = 4;
elseif P == 6
    Name = 'T1';
    x_L = 0;
    x_U = 10;
    x_0 = 0;
    x_opt = 5.22406;
    f_opt = -1.61642;
    f_Low=-500;
    x_min = 0;
    x_max = 10;
    nGlobal = 1;
    nLocal  = 5;
elseif P == 7
    Name = 'T2';
    x_L = -5;
    x_U = 2;
    x_0 = 0;
    x_opt = -2.84239;
    f_opt = -1.48911;
    f_Low=-500;
    x_min = -5;
    x_max = 2;
    nGlobal = 1;
    nLocal  = 5;
elseif P == 8
    Name = 'T3';
    x_L = -10;
    x_U = 10;
    x_0 = 0;
    x_opt = [-7.397286833804756; -1.114110069291653; 5.169069975025649];
    f_opt = -14.837950025122410; % f_opt = [-14.837950025122410; -14.837950006966516; -14.837949983145515]
    f_Low=-500;
    x_min = -10;
    x_max = 10;
    nGlobal = 1;
    nLocal  = 20;
elseif P == 9
    Name = 'T4';
    x_L = -5;
    x_U = 10;
    x_0 = 0;
    x_opt = 7.85411;
    f_opt = -0.999612;
    f_Low=-500;
    x_min = -5;
    x_max = 10;
    nGlobal = 1;
    nLocal  = 5;
elseif P == 10
    Name = 'T5';
    x_L = -10;
    x_U = 10;
    x_0 = 0;
    x_opt = -9.44999;
    f_opt = -0.943045;
    f_Low=-500;
    x_min = -10;
    x_max = 10;
    nGlobal = 1;
    nLocal  = 7;
elseif P == 11
    Name = 'T6';
    x_L = -5;
    x_U = 5;
    x_0 = 0;
    x_opt = 0;
    f_opt = -1;
    f_Low=-500;
    x_min = -5;
    x_max = 5;
    nGlobal = 1;
    nLocal  = 27;
elseif P == 12
    Name = 'T7';
    x_L = -30;
    x_U = 30;
    x_0 = 0;
    x_opt = 0;
    f_opt = 0;
    f_Low=-500;
    x_min = -30;
    x_max = 30;
    nGlobal = 1;
    nLocal  = 50;
elseif P == 13
    Name = 'T8';
    x_L = -1;
    x_U = 1;
    x_0 = 0;
    x_opt = 0;
    f_opt = -0.1;
    f_Low=-500;
    x_min = -1;
    x_max = 1;
    nGlobal = 1;
    nLocal  = 5;
elseif P == 14
    Name = 'T9';
    x_L = 0;
    x_U = 10;
    x_0 = 0;
    x_opt = 8.00754;
    f_opt = -.988082;
    f_Low=-500;
    x_min = 0;
    x_max = 10;
    nGlobal = 1;
    nLocal  = 7;
elseif P == 15
    Name = 'T10';
    x_L = -100;
    x_U = 100;
    x_0 = 0;
    x_opt = 0;
    f_opt = 0;
    f_Low=-500;
    x_min = -100;
    x_max = 100;
    nGlobal = 1;
    nLocal  = 30;
elseif P == 16
    Name = 'T11';
    x_L = -5;
    x_U = 5;
    x_0 = 0;
    x_opt = 1;
    f_opt = 0;
    f_Low=-500;
    x_min = -5;
    x_max = 5;
    nGlobal = 1;
    nLocal  = 25;
elseif P == 17
    Name = 'T12';
    x_L = 0;
    x_U = 10;
    x_0 = 0;
    x_opt = 9.57425;
    f_opt = -1.24825;
    f_Low=-500;
    x_min = 0;
    x_max = 10;
    nGlobal = 1;
    nLocal  = 3;
elseif P == 18
    Name = 'T13';
    x_L = -.5;
    x_U = 3;
    x_0 = 0;
    x_opt = 0.393797327485793 + (0:3)'; % Infinately many solutions (-N:N)'
    f_opt = -2.499607604320059;
    f_Low=-500;
    x_min = -.5;
    x_max = 3;
    nGlobal = 1;
    nLocal  = 10;
elseif P == 19
    Name = 'T14';
    x_L = 0;
    x_U = 20;
    x_0 = 0;
    x_opt =  7.853981154220392; % x_opt = [ 1.579547495907403;  7.853981154220392; 14.137177483360627];
    f_opt = -0.999999999941152; % f_opt = [-0.991134613568260; -0.999999999941152; -0.999999999833293];
    f_Low=-500;
    x_min = 0;
    x_max = 20;
    nGlobal = 1;
    nLocal  = 4;
elseif P == 20
    Name = 'T15';
    x_L = .2;
    x_U = 7;
    x_0 = 0;
    x_opt = 6.28319;
    f_opt = 15;
    f_Low=-500;
    x_min = -.2;
    x_max = 7;
    nGlobal = 1;
    nLocal  = 7;
elseif P == 21
    Name = 'T16';
    x_L = .2;
    x_U = 7;
    x_0 = 0;
    x_opt = 6.89453;
    f_opt = -17.5829;
    f_Low=-500;
    x_min = .2;
    x_max = 7;
    nGlobal = 1;
    nLocal  = 22;
elseif P == 22
    Name = 'T17';
    x_L = .2;
    x_U = 7;
    x_0 = 0;
    x_opt = 2.83935;
    f_opt = -0.952897;
    f_Low=-500;
    x_min = .2;
    x_max = 7;
    nGlobal = 1;
    nLocal  = 6;
elseif P == 23
    Name = 'T18';
    x_L = .2;
    x_U = 7;
    x_0 = 0;
    x_opt = 6.92006;
    f_opt = -6.26287;
    f_Low=-500;
    x_min = .2;
    x_max = 7;
    nGlobal = 1;
    nLocal  = 4;
elseif P == 24
    Name = 'T19';
    x_L = .2;
    x_U = 7;
    x_0 = 0;
    x_opt = 5.13434;
    f_opt = -7.04744;
    f_Low=-500;
    x_min = -.2;
    x_max = 7;
    nGlobal = 1;
    nLocal  = 3;
elseif P == 25
    Name = 'T20';
    x_L = .2;
    x_U = 7;
    x_0 = 0;
    x_opt = 4.85806;
    f_opt = 59.1291;
    f_Low = -500;
    x_min = .2;
    x_max = 7;
    nGlobal = 1;
    nLocal  = 8;
elseif P == 26
    Name = 'T21';
    x_L = -10;
    x_U = 30;
    x_0 = 0;
    x_opt = 12.855;
    f_opt = -1.73488;
    f_Low = -500;
    x_min = -10;
    x_max = 30;
    nGlobal = 1;
    nLocal  = 13;
elseif P == 27
    Name = 'T22';
    x_L = -5;
    x_U = 5;
    x_0 = 0;
    % Inifinately many global optima.
    x_opt = [-2.942197039745906; -0.199375584406395;  3.340969345255971]; 
    f_opt = -3.200092319187990;  % f_opt = [-3.200092312004187; -3.200092317694163; -3.200092319187990];
    f_Low = -500;
    x_min = -5;
    x_max = 5;
    nGlobal = 1;
    nLocal  = 17;
else
    error('lgo1_prob: Illegal problem number');
end

% Set x_0 to zeros
if isempty(x_0)
    x_0=zeros(length(x_L),1);
end

Prob = glcAssign('lgo1_f', x_L, x_U, Name, [], [], [], ...
    [], [], [], x_0, [], [], [], [], ...
    f_Low, x_min, x_max, f_opt, x_opt);
Prob.MIP.nLocal  = nLocal;
Prob.MIP.nGlobal = nGlobal;
Prob.P = P;

% MODIFICATION LOG:
%
% 040114  med  Created
% 041117  med  xxx_prob removed and code added
% 060928  med  Error updated
% 080603  med  Switched to glcAssign, cleaned
% 080916  nhq  Added multiple global optima for some problems