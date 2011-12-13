% Defines global optimization problems.
%
% lgo2_prob: Defines global optimization problems.
% Two or more dimensions
%
% The problems are defined by:
% Prob = probInit('lgo2_prob', i); i = 1,...,53
%
% function [probList, Prob] = lgo2_prob(P, ask, Prob);
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
% Written June 1, 2004.   Last modified Sept 29, 2008.

function [probList, Prob] = lgo2_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'M1'...
    ,'M2'...
    ,'M3'...
    ,'M4'...
    ,'M5'...
    ,'M6'...
    ,'M7'...
    ,'M8'...
    ,'M9'...
    ,'M10'...
    ,'M11'...
    ,'M12'...
    ,'M13'...
    ,'M14'...
    ,'M15'...
    ,'M16'...
    ,'M17'...
    ,'M18'...
    ,'M19'...
    ,'M20'...
    ,'M21'...
    ,'M22'...
    ,'M23'...
    ,'M24'...
    ,'M25'...
    ,'M26'...
    ,'M27'...
    ,'M28'...
    ,'M29'...
    ,'M30'...
    ,'M31'...
    ,'M32'...
    ,'M33'...
    ,'M34'...
    ,'M35'...
    ,'M36'...
    ,'M37'...
    ,'M38'...
    ,'M39'...
    ,'M40'...
    ,'M41'...
    ,'M42'...
    ,'M43'...
    ,'M44'...
    ,'M45'...
    ,'M46'...
    ,'M47'...
    ,'M48'...
    ,'M49'...
    ,'M50'...
    ,'M51'...
    ,'M52'...
    ,'M53'); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

if P == 1
    Name = 'M1';
    x_L = [0 0]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [5.22406 5.22406];
    f_opt = -3.23285;
    f_Low = -200;
    x_min = [0 0]';
    x_max = [10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 2
    Name = 'M2';
    x_L = [-5 -5]';
    x_U = [5 5]';
    x_0 = [0 0]';
    x_opt = [3.34098 -.199384];
    f_opt = -6.40018;
    f_Low = -500;
    x_min = [-5 -5]';
    x_max = [5 5]';
    nGlobal = 1;
    nLocal  = 2;
elseif P == 3
    Name = 'M3';
    x_L = [-30 -30]';
    x_U = [30 30]';
    x_0 = [0 0]';
    x_opt = [0 0];
    f_opt = 0;
    f_Low = -500;
    x_min = [-30 -30]';
    x_max = [30 30]';
    nGlobal = 1;
    nLocal  = 2;
elseif P == 4
    Name = 'M4';
    x_L = [-5 -5]';
    x_U = [5 5]';
    x_0 = [0 0]';
    x_opt = [-.691792 -2.82460];
    f_opt = -130.128;
    f_Low = -500;
    x_min = [-5 -5]';
    x_max = [5 5]';
    nGlobal = 1;
    nLocal  = 2;
elseif P == 5
    Name = 'M5';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [5 5];
    f_opt = 0;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 4;
    nLocal  = 4;
elseif P == 6
    Name = 'M6';
    x_L = [-5 -5]';
    x_U = [5 5]';
    x_0 = [0 0]';
    x_opt = [0 0];
    f_opt = 0;
    f_Low = -500;
    x_min = [-5 -5]';
    x_max = [5 5]';
    nGlobal = 1;
    nLocal  = 5;
elseif P == 7
    Name = 'M7';
    x_L = [-5 -5]';
    x_U = [5 5]';
    x_0 = [0 0]';
    x_opt = [0 0];
    f_opt = 0;
    f_Low = -500;
    x_min = [-5 -5]';
    x_max = [5 5]';
    nGlobal = 1;
    nLocal  = 5;
elseif P == 8
    Name = 'M8';
    x_L = [-5 -5]';
    x_U = [5 5]';
    x_0 = [0 0]';
    x_opt = [0 0];
    f_opt = 0;
    f_Low = -500;
    x_min = [-5 -5]';
    x_max = [5 5]';
    nGlobal = 1;
    nLocal  = 20;
elseif P == 9
    Name = 'M9';
    x_L = [-5 0]';
    x_U = [10 15]';
    x_0 = [0 0]';
    % x_opt = [-pi 12.275]; %NHQ, Multiple optima:
    x_opt = [ pi 2.275 ; 3*pi 2.475 ; -pi 12.275 ];
    f_opt = 5/(4*pi);
    f_Low = -500;
    x_min = [-5 0]';
    x_max = [10 15]';
    nGlobal = 1;
    nLocal  = 5;
elseif P == 10
    Name = 'M10';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [0.402537 0.287408];
    f_opt = 0;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 5;
    nLocal  = 5;
elseif P == 11
    Name = 'M11';
    x_L = [-5 -5]';
    x_U = [5 5]';
    x_0 = [0 0]';
    x_opt = [0 0];
    f_opt = 0;
    f_Low = -500;
    x_min = [-5 -5]';
    x_max = [5 5]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 12
    Name = 'M12';
    x_L = [-5 -5]';
    x_U = [5 5]';
    x_0 = [0 0]';
    x_opt = [0.089842 -0.712656 ; -0.089842 0.712656];
    f_opt = -1.031628453488552;
    f_Low = -500;
    x_min = [-5 -5]';
    x_max = [5 5]';
    nGlobal = 2;
    nLocal  = 2;
elseif P == 13
    Name = 'M13';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [1 1];
    f_opt = 0;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 1;
    nLocal  = 5;
elseif P == 14
    Name = 'M14';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [pi pi];
    f_opt = -1;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 15
    Name = 'M15';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [6.52448 -4.71239];     %NHQ, multiple optima??
    f_opt = -1.98816;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 16
    Name = 'M16';
    x_L = [-50 -50]';
    x_U = [50 50]';
    x_0 = [0 0]';
    x_opt = [0 0];
    f_opt = 0;
    f_Low = -500;
    x_min = [-50 -50]';
    x_max = [50 50]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 17
    Name = 'M17';
    x_L = [-2 -2]';
    x_U = [2 2]';
    x_0 = [0 0]';
    x_opt = [0 -1];
    f_opt = 3;
    f_Low = -500;
    x_min = [-2 -2]';
    x_max = [2 2]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 18
    Name = 'M18';
    x_L = [-100 -100]';
    x_U = [100 100]';
    x_0 = [0 0]';
    x_opt = [0 0];
    f_opt = 0;
    f_Low = -500;
    x_min = [-100 -100]';
    x_max = [100 100]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 19
    Name = 'M19';
    x_L = [-100 -100]';
    x_U = [100 100]';
    x_0 = [0 0]';
    x_opt = [5.93234 -1.30671];     %NHQ, multiple optima??
    f_opt = -137.5396;
    f_Low = -500;
    x_min = [-100 -100]';
    x_max = [100 100]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 20
    Name = 'M20';
    x_L = [0 0]';
    x_U = [5 5]';
    x_0 = [0 0]';
    x_opt = [4 2];
    f_opt = -2.34581;
    f_Low = -500;
    x_min = [0 0]';
    x_max = [5 5]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 21
    Name = 'M21';
    x_L = [-5 -5]';
    x_U = [5 5]';
    x_0 = [0 0]';
    x_opt = [1 1];
    f_opt = 0;
    f_Low = -500;
    x_min = [-5 -5]';
    x_max = [5 5]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 22
    Name = 'M22';
    x_L = [-5 -5]';
    x_U = [5 5]';
    x_0 = [0 0]';
    x_opt = [1 1];
    f_opt = 0;
    f_Low = -500;
    x_min = [-5 -5]';
    x_max = [5 5]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 23
    Name = 'M23';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [-1 -1];
    f_opt = 0;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 24
    Name = 'M24';
    x_L = [-3 4.1]';
    x_U = [12.1 5.8]';
    x_0 = [0 0]';
    x_opt = [11.6255446616 5.7250442438];
    f_opt = -38.8502944794;
    f_Low = -500;
    x_min = [-3 4.1]';
    x_max = [12.1 5.8]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 25
    Name = 'M25';
    x_L = [-2 -2]';
    x_U = [2 2]';
    x_0 = [0 0]';
    x_opt = [-0.0135407 -0.0135407];
    f_opt = -1.29695;
    f_Low = -500;
    x_min = [-2 -2]';
    x_max = [2 2]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 26
    Name = 'M26';
    x_L = [-4 -4]';
    x_U = [4 4]';
    x_0 = [0 0]';
    x_opt = [2 2];
    f_opt = -2;
    f_Low = -500;
    x_min = [-4 -4]';
    x_max = [4 4]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 27
    Name = 'M27';
    x_L = [-15 -15]';
    x_U = [15 15]';
    x_0 = [0 0]';
    x_opt = [1.041268 1.341268];
    f_opt = -1.0459;
    f_Low = -500;
    x_min = [-15 -15]';
    x_max = [15 15]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 28
    Name = 'M28';
    x_L = [3 3]';
    x_U = [9.999 9.999]';
    x_0 = [0 0]';
    x_opt = [8.53879 8.53879];
    f_opt = 4.98151;
    f_Low = -500;
    x_min = [3 3]';
    x_max = [9.999 9.999]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 29
    Name = 'M29';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [0 0];
    f_opt = 0;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 30
    Name = 'M30';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [1 1];
    f_opt = 0;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 31
    Name = 'M31';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [1.14030 0.89841];
    f_opt = 1.95259;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 32
    Name = 'M32';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [1 1];
    f_opt = 2;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 33
    Name = 'M33';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [1.2 2.4];
    f_opt = 7.2;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 34
    Name = 'M34';
    x_L = [-10 -10 -10 -10]';
    x_U = [10 10 10 10]';
    x_0 = [0 0 0 0]';
    x_opt = [0 1 2 -1];
    f_opt = -44;
    f_Low = -500;
    x_min = [-10 -10 -10 -10]';
    x_max = [10 10 10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 35
    Name = 'M35';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [0.5 0.5];
    f_opt = 0;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 36
    Name = 'M36';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [1 1];
    f_opt = 0;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 37
    Name = 'M37';
    x_L = [-10 -10 -10 -10]';
    x_U = [10 10 10 10]';
    x_0 = [0 0 0 0]';
    x_opt = [1 1 1 1];
    f_opt = 0;
    f_Low = -500;
    x_min = [-10 -10 -10 -10]';
    x_max = [10 10 10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 38
    Name = 'M38';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    %NHQ, this is a saddle point, but not the global min.
%     x_opt = [0.5 0.5];
%     f_opt = 0;
    x_opt = [1.025 -10];
    f_opt = -7.35;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 39
    Name = 'M39';
    x_L = [-3 -3]';
    x_U = [3 3]';
    x_0 = [0 0]';
    x_opt = [0 0];
    f_opt = -2;
    f_Low = -500;
    x_min = [-3 -3]';
    x_max = [3 3]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 40
    Name = 'M40';
    x_L = [-5 -5]';
    x_U = [5 5]';
    x_0 = [0 0]';
    x_opt = [0 0];
    f_opt = 0;
    f_Low = -500;
    x_min = [-5 -5]';
    x_max = [5 5]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 41
    Name = 'M41';
    x_L = [0.5 0.5]';
    x_U = [5 5]';
    x_0 = [0 0]';
    x_opt = [0.5 1.5];
    f_opt = 0;
    f_Low = -500;
    x_min = [0.5 0.5]';
    x_max = [5 5]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 42
    Name = 'M42';
    x_L = [1 1]';
    x_U = [5 5]';
    x_0 = [0 0]';
    %x_opt = [pi/2 pi];      %NHQ, multiple optima??
    x_opt = [ ];
    for i = 1:3
       for j = 1:4
          x_opt = [x_opt; [i*pi/2 j*pi/3]];
       end
    end
    f_opt = 0;
    f_Low = -500;
    x_min = [1 1]';
    x_max = [5 5]';
    nGlobal = 2;
    nLocal  = 2;
elseif P == 43
    Name = 'M43';
    x_L = [-100 -100]';
    x_U = [100 100]';
    x_0 = [0 0]';
    x_opt = [0 0];
    f_opt = 0;
    f_Low = -500;
    x_min = [-100 -100]';
    x_max = [100 100]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 44
    Name = 'M44';
    x_L = [-100 -100]';
    x_U = [100 100]';
    x_0 = [0 0]';
    x_opt = [0 0];
    f_opt = 0;
    f_Low = -500;
    x_min = [-100 -100]';
    x_max = [100 100]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 45
    Name = 'M45';
    x_L = [-100 -100]';
    x_U = [100 100]';
    x_0 = [0 0]';
    x_opt = [0 0];
    f_opt = 0;
    f_Low = -500;
    x_min = [-100 -100]';
    x_max = [100 100]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 46
    Name = 'M46';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [-7.0835 4.8580];       %NHQ, multiple optima??
    f_opt = -186.7309;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 18;
    nLocal  = 760;
elseif P == 47
    Name = 'M47';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [-1.42513 -0.80032];    %NHQ, multiple optima??
    f_opt = -186.730909;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 48
    Name = 'M48';
    x_L = [-10 -10]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [-6.774576 -0.491391];  %NHQ, multiple optima??
    f_opt = -24.062499;
    f_Low = -500;
    x_min = [-10 -10]';
    x_max = [10 10]';
    nGlobal = 9;
    nLocal  = 9;
elseif P == 49
    Name = 'M49';
    x_L = [-500 -500]';
    x_U = [500 500]';
    x_0 = [0 0]';
    x_opt = [420.97 420.97];
    f_opt = -837.9658;
    f_Low = -500;
    x_min = [-500 -500]';
    x_max = [500 500]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 50
    Name = 'M50';
    x_L = [0 0]';
    x_U = [180 180]';
    x_0 = [0 0]';
    x_opt = [113.252205 78.694686];
    f_opt = -3.5;
    f_Low = -500;
    x_min = [0 0]';
    x_max = [180 180]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 51
    Name = 'M51';
    x_L = [0 0]';
    x_U = [10 10]';
    x_0 = [0 0]';
    x_opt = [8.02406531 9.14653404];
    f_opt = -12.11900837975;
    f_Low = -500;
    x_min = [0 0]';
    x_max = [10 10]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 52
    Name = 'M52';
    x_L = [-3 -3]';
    x_U = [3 3]';
    x_0 = [0 0]';
    x_opt = [-0.0244030796 0.2106124272];
    f_opt = -3.306868647;
    f_Low = -500;
    x_min = [-3 -3]';
    x_max = [3 3]';
    nGlobal = 1;
    nLocal  = 1;
elseif P == 53
    Name = 'M53';
    x_L = [-2 -2]';
    x_U = [2 2]';
    x_0 = [0 0]';
    x_opt = [0 0];
    f_opt = -2;
    f_Low = -500;
    x_min = [-2 -2]';
    x_max = [2 2]';
    nGlobal = 1;
    nLocal  = 1;
else
    disp('lgo2_prob: Illegal problem number')
    pause
    Name=[];
end
% Set x_0 to zeros (dummy for GUI)
if isempty(x_0)
    x_0=zeros(length(x_L),1);
end

Prob = glcAssign('lgo2_f', x_L, x_U, Name, [], [], [], ...
    [], [], [], x_0, [], [], [], [], ...
    f_Low, x_min, x_max, f_opt, x_opt);
Prob.MIP.nLocal  = nLocal;
Prob.MIP.nGlobal = nGlobal;
Prob.P = P;

% MODIFICATION LOG:
%
% 040114  med  Created
% 040204  hkh  Corrected wrong solution for problem 27
% 041117  med  xxx_prob removed and code added
% 060928  med  Error updated
% 080603  med  Switched to glcAssign, cleaned'
% 080918  nhq  Added alternative optimas for some problems.
% 080929  hkh  M16, M17 missing in list of names
