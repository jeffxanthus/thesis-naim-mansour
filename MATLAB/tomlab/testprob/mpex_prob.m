% mpexp_prob:
%
% Defines Mixed Complementarity Programming example problems
%
% function [probList, Prob] = mpex_prob(P);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Anders Goran, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2006-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written Jun 22, 2006.  Last modified Jun 3, 2008.

function [probList, Prob] = mpex_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'Gauvin-Savard' ...
    ,'DF1 - QQR2-MN-4-2' ...
    ,'BiLin QQR2-MN-8-5' ...
    ,'Bard3m' ...
    ,'RalphMOD'...
    ); % CHANGE: MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return;
end

HessPattern = []; ConsPattern = []; c_L = []; c_U = [];
pSepFunc = []; f_Low = []; x_opt = []; f_opt = [];
A = []; x_min = []; x_max = []; b_L = []; b_U = []; x_0 = [];
user = [];

if P == 1
    % /* -------------------------------------------------------
    %
    %   MPEC example taken from Gauvin and Savard, "Working Paper G9037",
    %   GERAD, Ecole Polytechnique de Montreal (1992) (first version).
    %
    %    Let Y(x) := { y | 0 <= y <= 20 - x }.
    %    min  f(x,y) := x**2 + (y-10)**2
    %    s.t.      0 <= x <= 15
    %              y solves MCP((F(y) := 4(x + 2y - 30)), Y(x))
    %
    %   We have modified the problem by adding a dual variable and
    %   incorporating the constraint y <= 20 - x into the MCP function F.
    %
    %   From a GAMS model by S.P. Dirkse & M.C. Ferris (MPECLIB),
    %   (see http://www.gams.com/mpec/).
    %
    %   AMPL coding Sven Leyffer, University of Dundee, Jan. 2000
    %
    %   ------------------------------------------------------- */
    %
    % var x >= 0, <= 15; # design variable
    % var y >= 0;		 # state variable
    % var u >= 0;        # duals in MCP to determine state vars
    x_L = [0,0,0]';
    x_U = [15,inf,inf]';
    % minimize theta: x^2 + (y-10)^2;
    %
    % subject to
    %
    %    # ... F(x,y) from original problem
    %    Fy: 0 <= 4 * (x + 2*y - 30) + u   complements   y >= 0;
    %    Fu: 0 <=  20 - x - y              complements   u >= 0;
    Name = 'Gauvin-Savard';
    x_0 = [7.5, 1, 0]';
    c_L = [0,0]';
    c_U = [inf,inf];
    ConsPattern = [];
    HessPattern = [];
    % Complementarity: x(2) _|_ c(1)
    % and              x(3) _|_ c(2)
    mpec = [ ...
        2,0, 0,0, 1,0 ; ...
        3,0, 0,0, 2,0 ; ];
elseif P == 2
    % An MPEC from S.P. Dirkse and M.C. Ferris, Modeling & Solution
    % Environment for MPEC: GAMS & MATLAB, University of Wisconsin
    % CS report, 1997.
    %
    % min ( x(1) - 1 - x(2) )^2
    %
    % s/t
    %       x(1)^2                 <= 2
    %      (x(1)-1)^2 + (x(2)-1)^2 <= 3
    %
    %   0 <= x(2) - x(1)^2 + 1   complements x(2) >= 0;
    Name = 'DF1 - QQR2-MN-4-2';
    x_L = [-1; 0];
    x_U = [ 2, inf];
    % Third constraint is the first component in the complementarity.
    c_L = [-inf; -inf ; 0   ];
    c_U = [ 2  ; 3    ; inf ];
    % variable two versus nonlinear constraint three:
    mpec = [ 2,0, 0,0, 3,0 ];
elseif P == 3
    % An bilevel linear program due to Hansen, Jaumard and Savard,
    % "New branch-and-bound rules for linear bilevel programming,
    % SIAM J. Sci. Stat. Comp. 13, 1194-1217, 1992. See also book
    % Mathematical Programs with Equilibrium Constraints,
    % by Luo, Pang & Ralph, CUP, 1997, p. 357.
    %
    % maximize f: 8*x(1) + 4*x(2) - 4*x(3) + 40*x(4) + 4*x(5);
    %
    % subject to
    %
    %    lin:    x(1) + 2*x(2) - x(5) <= 1.3;
    %
    %    KKT1:   0 <= 2 - x(6) - 2*x(7) + 4*x(8)           complements  x(3) >= 0;
    %    KKT2:   0 <= 1 + x(6) + 4*x(7) - 2*x(8)           complements  x(4) >= 0;
    %    KKT3:   0 <= 2 + x(6) - x(7) - x(8)               complements  x(5) >= 0;
    %
    %    slack1: 0 <= 1 + x(3) - x(4) - x(5)               complements  x(6) >= 0;
    %    slack2: 0 <= 2 - 4*x(1) + 2*x(3) - 4*x(4) + x(5)  complements  x(7) >= 0;
    %    slack3: 0 <= 2 - 4*x(2) - 4*x(3) + 2*x(4) + x(5)  complements  x(8) >= 0;
    mpec = [ ...
        3,0,2,0,0,0, ; ...
        4,0,3,0,0,0, ; ...
        5,0,4,0,0,0, ; ...
        6,0,5,0,0,0, ; ...
        7,0,6,0,0,0, ; ...
        8,0,7,0,0,0, ; ...
        ];
    x_L = zeros(8,1);
    x_U = inf*ones(8,1);
    x_0 = 1.0*ones(8,1);
    A = [ ...
        1     2     0     0    -1     0     0     0   ; ...
        0     0     0     0     0    -1    -2     4   ; ...
        0     0     0     0     0     1     4    -2   ; ...
        0     0     0     0     0     1    -1    -1   ; ...
        0     0     1    -1    -1     0     0     0   ; ...
        -4     0     2    -4     1     0     0     0   ; ...
        0    -4    -4     2     1     0     0     0   ; ...
        ];
    b_L = [-inf,  -2,  -1,  -2,  -1,  -2,  -2 ];
    b_U = [ 1.3, inf, inf, inf, inf, inf, inf ];
    Name = 'bilin QQR2-MN-8-5';
elseif P == 4
    Name = 'Bard3m';
    x_L = zeros(6,1);
    x_U = inf*ones(6,1);
    % # .. upper level problem objective
    % minimize cost: -x1^2 - 3*x2 + x4^2 - 4*x3;
    %
    % subject to
    %    side: x1^2 + 2*x2 <= 4;
    %
    %    # ... optimality conditions of lower level problem
    %    cons1: 0 <= x1^2 - 2*x1 + x2^2 - 2*x3 + x4 + 3   complements  x5 >= 0;
    %    cons2: 0 <= x2 + 3*x3 - 4*x4 - 4                 complements  x6 >= 0;
    %
    %    d_x3 : 0 <= (2*x3+2*x5)-3*x6   complements  x3 >= 0;
    %    d_x4 : 0 <= (-5-x5)+4*x6       complements  x4 >= 0;
    c_L = [-inf, 0,0,0,0 ];
    c_U = [ 4, inf, inf, inf, inf ];
    % Complementarities:
    mpec = [ ...
        5,0, 0,0, 2,0 ; ...
        6,0, 0,0, 3,0 ; ...
        3,0, 0,0, 4,0 ; ...
        4,0, 0,0, 5,0 ];
    f_opt = -12.6787109375;
elseif P == 5
    Name = 'Ralph';
    % Ralphmod.mat contains a different variable P...
    Psave = P;
    load Ralphmod.mat M P q c
    user.P = P;
    user.c = c;
    P = Psave;
    nd = 4;    % 4 design variables
    ns = 100;  % 100 state variables
    N = nd+ns;
    A   = M(nd+1:end, :);
    b_L = -q(nd+1:end);
    b_U = inf*ones(size(b_L));
    x_L           = zeros( N, 1);
    x_U           = ones( N, 1);
    x_U(1:nd)     = 1000;  % Design variables <= 1000
    x_U(nd+1:end) = inf;   % State variables unlimited above
    c_L = [];
    c_U = [];
    % The last ns variables are complementary to the (first) ns linear
    % cosntraints:
    mpec = [(nd+1:N)', zeros(ns,1), (1:ns)', zeros(ns,3)];
    f_opt = -683.033;
else
    error('mpex_prob: Illegal problem number');
end

% Define the Prob
c  = 'mpex_c';
dc = 'mpex_dc';

if isempty(c_L) & isempty(c_U)
    c  = [];
    dc = [];
end

% Define the Prob
Prob = conAssign('mpex_f','mpex_g',[], HessPattern, x_L,...
    x_U, Name, x_0, pSepFunc, f_Low, A, b_L, b_U, c, dc,...
    [], ConsPattern, c_L, c_U, x_min, x_max, f_opt, x_opt);
Prob.P = P;
Prob.user = user;
Prob   = BuildMPEC(Prob,mpec);

% MODIFICATION LOG:
%
% 060622 ango Wrote file
% 080603 med  Switched to conAssign, cleaned