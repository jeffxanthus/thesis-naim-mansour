% minlp_prob.m
%
% Defines mixed-integer nonlinear programming (MINLP) problems.
%
% function [probList, Prob] = minlp_prob(P);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2009 by Tomlab Optimization Inc. $Release: 7.3.0$
% Written June 1, 1999.   Last modified Oct 8, 2009.

function [probList, Prob] = minlp_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'Kocis & Grossmann 1988'...
    ,'Floudas 1995 6.6.5'...
    ,'Yuan et al. 1988'...
    ,'Poern et al. 1997'...
    ,'Trim loss 1'...
    ,'Trim loss 2'...
    ,'Trim loss 3'...
    ,'Trim loss 4'...
    ,'Yuan et al. 1988 (Alt)'...
    ,'Kocis & Grossmann 1988 (Alt)'...
    ,'Kocis & Grossmann 1989'...
    ,'Lee & Grossmann Disjunctive B&B'...
    ,'Kesavan et al. 2004 B'...
    ,'Kesavan et al. 2004 D'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

HessPattern = []; VarWeight = []; ConsPattern = [];
user = [];
% If to run CGO tests
BIGBND = 1E4;
BIG    = inf;

if P == 1
    Name = 'Kocis & Grossmann 1998';
    IntVars = [3 4 5];
    % May get division by zero if not x_L(2)==small but >0
    x_L = [ 0      1/BIGBND   0 0 0 ]';
    x_U = [ BIGBND  BIGBND    1 1 1 ]';
    A = [1 0     1  0 0 ; ...
        0 1.333 0  1 0 ; ...
        0 0    -1 -1 1 ];
    b_L = -inf*ones(3,1);
    b_U = [1.6 ; 3 ; 0];
    c_L = [1.25;3];
    c_U = c_L;
    % Setting x_0 = ones(5,1); gives a non-feasible point with respect to
    % linear constraints. Adjusting this point gives (0.6,1,1,1,1)
    % However from this initial point MINLP will converge to a
    % nonglobal solution
    %
    % If disabling the feasibility test, then MINLP converges
    x_0 = ones(5,1);
    x_opt = [1.12,1.31,0,1,1];
    f_opt = 7.6672;
    x_min = [-1 -1 0 0 0];
    x_max = [ 1  1 1 1 1];
    ConsPattern = [1 0 1 0 0 ; 0 1 0 1 0];
elseif P == 10
    % This is the same problem as problem 1, but the linear constraints
    % are put among the nonlinear.
    % Thereby MINLPBB never checks for feasibility of the linear constraints
    % and converge to the global solution from the initial point x=(1,1,1,1,1)
    Name = 'Kocis & Grossmann 1998-X';
    IntVars = [3 4 5];
    % May get division by zero if not x_L(2)==small but >0
    x_L = [ 0      1/BIGBND   0 0 0 ]';
    x_U = [ BIGBND  BIGBND    1 1 1 ]';
    b_L = [-inf;-inf;-inf];
    b_U = [1.6 ; 3 ; 0];
    c_L = [1.25;3];
    c_U = c_L;
    c_L = [c_L;b_L];
    c_U = [c_U;b_U];
    x_0 = ones(5,1);
    x_opt = [1.118,1.310,0,1,1];
    f_opt = 7.66718;
    x_min = [-1 -1 0 0 0];
    x_max = [ 1  1 1 1 1];
    A   = [];
    b_L = [];
    b_U = [];
elseif P==2
    Name = 'Floudas 1995 6.6.5';
    % glc_prob: P=20: Floudas-Pardalos 12.2 TP 2
    x_L = [0.2; -2.22554 ; 0  ];
    x_U = [1.0; -1.0     ; 1  ];
    IntVars = logical([0 0 1]);
    A = [ 0.0  1.0  1.1 ; ...
        1.0  0.0 -1.2 ];
    b_L = [ -inf ; -inf ];
    b_U = [  -1  ;  0.2 ];
    c_L = -inf;
    c_U = 0;
    f_opt = 1.076543;
    x_opt = [ 0.941937 ; -2.1 ; 1 ];
    x_0 = x_opt;
    x_min = x_L;
    x_max = x_U;
elseif P==3
    Name = 'Yuan et al. 1988';
    % glc_prob: P == 21, Name = 'Floudas-Pardalos 12.2 TP 3';
    IntVars = 1:4;
    x_L = zeros(7,1);
    x_U = [1 1 1 1 inf inf inf]';
    x_0 = x_L;
    A = [ 1 1 1 0 1 1 1 ; ...
        1 0 0 0 1 0 0 ; ...
        0 1 0 0 0 1 0 ; ...
        0 0 1 0 0 0 1 ; ...
        0 0 0 1 1 0 0 ];
    b_L = -BIG*ones(5,1);
    b_U = [5.0 ; 1.2 ; 1.8 ; 2.5 ; 1.2];
    c_L = -BIG*ones(4,1);
    c_U = [ 5.5 ; 1.64 ; 4.25 ; 4.64 ];
    x_opt = [1 1 0 1 .2 .8 1.907879];
    f_opt = 4.579582;
    ConsPattern = [ ...
        0 0 1 0 1 1 1 ; ...
        0 1 0 0 0 1 0 ; ...
        0 0 1 0 0 0 1 ; ...
        0 1 0 0 0 0 1 ];
    x_min = [0 0 0 0 0 0 0 ];
    x_max = [1 1 1 1 1 1 1 ];
elseif P==9
    % This is the same problem as problem 3, but the linear constraints
    % are put among the nonlinear.
    Name = 'Yuan et al. 1988-X';
    IntVars = 1:4;
    x_L = zeros(7,1);
    x_U = [1 1 1 1 inf inf inf]';
    x_0 = x_L;
    A = [ 1 1 1 0 1 1 1 ; ...
        1 0 0 0 1 0 0 ; ...
        0 1 0 0 0 1 0 ; ...
        0 0 1 0 0 0 1 ; ...
        0 0 0 1 1 0 0 ];
    b_L = -BIG*ones(5,1);
    b_U = [5.0 ; 1.2 ; 1.8 ; 2.5 ; 1.2];
    c_L = -BIG*ones(4,1);
    c_U = [ 5.5 ; 1.64 ; 4.25 ; 4.64 ];
    c_L = [c_L;b_L];
    c_U = [c_U;b_U];
    x_opt = [1 1 0 1 .2 .8 1.907879];
    f_opt = 4.5795815;
    ConsPattern = [ ...
        0 0 1 0 1 1 1 ; ...
        0 1 0 0 0 1 0 ; ...
        0 0 1 0 0 0 1 ; ...
        0 1 0 0 0 0 1;A ];
    A   = [];
    b_L = [];
    b_U = [];
    x_min = [];
    x_max = [];
elseif P==4
    Name = 'Pörn et al. 1997';
    % glc_prob, P == 23, Name = 'Floudas-Pardalos 12.2 TP 5';
    IntVars = 1:2;
    x_L = [1 1]';
    x_U = [5 5]';
    x_0 = x_L;
    A = [ ...
        -1 -2 ; ...
        -3  1 ; ...
         4 -3 ];
    % Bounds as given in Chapter 12, Handbook of Test Problems.
    % Bounds as given originally - gives different solution [1 1];
    b_L = -BIG*ones(3,1);
    b_U = [-5 1 11 ]';
    c_U = -24;
    c_L = -inf;
    ConsPattern = [];
    f_opt = 31;
    x_opt = [3 1];
    x_min = [];
    x_max = [];
elseif P>=5 & P<=8
    % The trim loss problems 1-4
    x_min = [];
    x_max = [];
    switch(P)
        case 5,
            N = 4;
            p = 4;
            Nkmax = 5;
            Bmax  = 1850;
            delta = 100;
            n = [ 15   28  21  30 ]';
            b = [ 290 315 350 455 ]';
            M = [ 30  30   30  30 ]';
            x_opt = [9 7 3 0  1 1 1 0  1 0 2 0 2 1 1 0 0 3 0 0 2 1 2 0];
            f_opt = 19.6;
        case 6,
            N = 4;
            p = 4;
            Nkmax = 5;
            Bmax  = 1900;
            delta = 200;
            n = [9 7 12 11]';
            b = [330 360 385 415]';
            M = [15 12 9 6]';
            % This solution from book gives f_opt = 11.1, obvious errors in book
            x_opt = [11 0 0 0  1 0 0 0  1 0 0 0  1 0 0 0  2 0 0 0  1 0 0 0];
            % The following solution gives 8.6, but is not feasible
            % x_opt = [8 0 0 0  1 1 1 0  1 0 0 0  1 0 0 0  2 0 0 0  1 0 0 0];
            %
            % minlpSolve found optimal feasible solution at 8.6 after 4960 iter
            % The optimal solution is:
            x_opt = [4 3 1 0  1 1 1 0  1 2 0 0  0 1 4 0  3 0 0 0  1 2 1 0];
            f_opt = 8.6;
        case 7,
            N = 5;
            p = 5;
            Nkmax = 5;
            Bmax  = 2000;
            delta = 200;
            n = [12 6 15 6 9]';
            b = [330 360 370 415 435]';
            M = [15 12 9 6 6]';
            % This solution from book gives f_opt = 15.1, obvious errors in book
            % Only length 30, should be 35
            % x_opt = [15 0 0 0 0  1 0 0 0 0  ...
            %     1 0 0 0 0  1 0 0 0 0  1 0 0 0 0  1 0 0 0 0];
            %
            % MINLPBB:
            % From "optimal" point
            % Sparse MINLP finds 11.5, optimal
            % Dense  MINLP says root infeasible
            % From x_0
            % Dense  MINLP finds 14.5, optimal
            % Sparse MINLP finds 11.5, optimal
            %
            % minlpSolve finds 10.3 after 977 iterations
            % The optimal solution is:
            x_opt = [6 4 0 0 0  1 1 0 0 0 ...
                  2 0 0 0 0  1 0 0 0 0  0 4 0 0 0  1 0 0 0 0  1 1 0 0 0]
            f_opt = 10.3;
        case 8,
            N = 6;
            p = 6;
            Nkmax = 5;
            Bmax  = 2200;
            delta = 100;
            n = [8 16 15 7 14 16]';
            b = [330 360 380 430 490 530]';
            M = [15 12 8 7 4 2]';
            % x_opt from book, f=15.3. One linear constraint is infeasible
            % Also nonlinear constraint 3 and 6 infeasible
            % f_opt = 15.3;
            % constraint 37: b_L=16    A*x_opt=15   b_U=Inf
            %x_opt = [8 7 0 0 0 0  1 1 0 0 0 0  1 0 0 0 0 0  2 0 0 0 0 0  ...
            %         0 2 0 0 0 0  0 1 0 0 0 0  0 2 0 0 0 0  1 0 0 0 0 0];
            %
            % minlpSolve found optimal solution 16.3 after 356 iter
            x_opt = [8 8 0 0 0 0  1 1 0 0 0 0  1 0 0 0 0 0  2 0 0 0 0 0  ...
                      0 2 0 0 0 0  0 1 0 0 0 0  0 2 0 0 0 0  2 0 0 0 0 0];
            f_opt = 16.3;
    end % switch(P)

    % Common for all four Trim Loss problems
    Name = ['Trim Loss ' num2str(P)-4 ];
    IntVars = (1:2*p + N*p);
    x_L = zeros( 2*p+N*p,1 );
    x_U = [M(:);ones(p,1);Nkmax*ones(N*p,1)];
    [A,b_L,b_U,c_L,c_U]=setup_trim(N,p,Nkmax,Bmax,delta,n,b,M,Inf);
    n = length(x_L);
    r = reshape(ones(n-2*p,1),N,p);
    ConsPattern = sparse([ r zeros(N,p) ]);
    for k=1:N
        ConsPattern(k, (k-1)*p+(2*p+1:3*p))=ones(1,p);
    end
    A = sparse(A);
    x_0 = [];
    % The following x_0 makes the problems converge for 5 and 8, sparse
    if P==5
        % minlpBB
        % This starting value gives integer infeasible, and f=21.3 for DENSE
        % And sometimes f=20.3 ... weird
        % Sparse converges fine ...
        %
        % minlpSolve converges fine to 19.6 in 2565 iterations
        x_0 = ones(length(x_L),1);
        x_0(1:length(M)) = round(0.5*M);
    elseif P==6
        % Book says f=8.6, but given optimal x has f(x)=11.1
        % MINLPBB runs:
        % If starting with x_opt with f=11.1, both dense and sparse stops.
        % Sparse finds f=8.6
        % Dense finds f=8.6, but says it is integer infeasible
        x_0 = 0.5*(x_U+x_L);
        % Dense finds optimal f=8.6
        % Sparse finds f=8.6, but claims stack overflow ...
        %
        % minlpSolve finds optimal solution f=8.6 in 4960 iterations
        x_0(1:length(M)) = M;
        x_0(1:length(M)) = round(0.5*M);
    elseif P==7
        % Book says f=10.3, but gives solution vector with f=15.1
        % Starting with given solution vector, sparse MINLPBB stops with 15.1
        % Starting with given solution vector, dense  MINLPBB stops with 15.1
        % BUT - sometimes with other x_0 sparse and dense finds f=10.6 as well
        % After presolve, dense finds 10.3, sparse finds 10.6
        %
        % minlpSolve finds 10.3 after 977 iterations
        x_0 = zeros(length(x_L),1);
        x_0(1:length(M)) = M;
    elseif P==8
        % Book says f=15.3, optimal solution is 16.3 (minlpSolve)
        % x_opt is linear and nonlinear infeasible
        % Runs with MINLPBB:
        % If start with x_opt, sparse and dense says stack overflow, f=16.3
        % Sparse gives stack overflow, finds 16.3
        % Dense says optimal solution found, but finds 16.3
        % Dense and sparse gives stack overflow, finds 16.3
        %
        % minlpSolve found optimal solution 16.3 after 356 iter
        x_0 = x_L;
        x_0(end-5:end) = 1;
        x_0(1:length(M)) = round(0.7*M);
    end
    user.N = N;
    user.P = p;
elseif P==11
    Name = 'Kocis & Grossmann 1989';
    % Variables [v,x,c,p,y1,y2];
    IntVars = [3 4];
    x_L = [  0  0  0 0 ]';
    x_U = [ 10 20  1 1 ]';
    A   = [0 0 1 1 ];
    b_L = 1;
    b_U = 1;
    c_L = -inf*ones(4,1);
    user.M = 100; % use big-M method
    x_0 = [0 0 0 1];
    % minlpSolve optimal solution:
    x_opt = [3.514243,13.427982,1,0];
    f_opt = 99.2396083;
    c_U   = [0;0;1;0];
    x_min = [ 0  0 0 0 ];
    x_max = [10 20 1 1 ];
    ConsPattern = [1 1 1 0 ; ...
        1 1 1 0 ; ...
        1 1 0 1 ; ...
        1 1 0 1];
elseif P==12
    Name = 'Lee & Grossmann Disjunctive B&B';
    % Variables [x(1),x(2),y1,y2,y3];
    IntVars = [3 4 5];
    x_L = [  0  0  0 0 0 ]';
    x_U = [  8  8  1 1 1 ]';
    A   = [0 0 1 1 1 ];
    b_L = 1;
    b_U = 1;
    c_L = -inf*ones(3,1);
    c_U = [0;0;0];
    user.M = 30; % use big-M method
    x_0 = [0 0 0 0 1];
    % minlpSolve optimal solution:
    x_opt = [3.292896,1.707112, 0, 1, 0];
    f_opt = 1.171572;
    x_min = x_L;
    x_max = x_U;
    ConsPattern = [1 1 1 0 0 ; ...
        1 1 0 1 0 ; ...
        1 1 0 0 1];
elseif P==13
    Name = 'Kesavan et al. 2004 B';
    IntVars = (1:8)';
    x_L = [0 0 0 0 0 0 0 0      0      0      0]';
    x_U = [1 1 1 1 1 1 1 1 0.9970 0.9985 0.9988]';
    A   = [-1 -1 -1  0  0  0  0  0  0  0  0;...
        0  0  0 -1 -1 -1  0  0  0  0  0;...
        0  0  0  0  0  0 -1 -1  0  0  0;...
        3  1  2  3  2  1  3  2  0  0  0];
    b_L = -inf*ones(4,1);
    b_U = [-1 -1 -1 10]';
    c_L = [0;0;0];
    c_U = [0;0;0];
    x_0 = [];
    % minlpSolve optimal solution:
    x_opt = [0 1 1 1 0 1 1 0 0.97 0.9925 0.98];
    f_opt = -0.9434706;
    x_min = x_L;
    x_max = x_U;
    ConsPattern = [1 1 1 0 0 0 0 0 1 0 0; ...
        0 0 0 1 1 1 0 0 0 1 0; ...
        0 0 0 0 0 0 1 1 0 0 1];
elseif P==14
    Name = 'Kesavan et al. 2004 D';
    IntVars = (1:3)';
    x_L = [0 0 0  1 1]';
    x_U = [1 1 1 10 6]';
    A   = [0 0 0 1 -1;...
        0 0 0 3  2;...
        1 2 4 0 -1];
    b_L = [-inf; -inf; 0];
    b_U = [3; 24; 0];
    c_L = -inf*ones(1,1);
    c_U = 39;
    x_0 = [0 0 0 1 1];
    % minlpSolve optimal solution:
    x_opt = [1 0 0 4 1];
    f_opt = -17;
    x_min = x_L;
    x_max = x_U;
    ConsPattern = [0 0 0 1 1];
else
    error('minlp_prob: Illegal problem number');
end

% Set x_0 to zeros (dummy for GUI)
if isempty(x_0)
    x_0=zeros(length(x_L),1);
end

Prob = minlpAssign('minlp_f','minlp_g','minlp_H', HessPattern, x_L, x_U, Name, x_0, ...
    IntVars, VarWeight, [], [], A, b_L, b_U, ...
    'minlp_c','minlp_dc','minlp_d2c', ConsPattern, c_L, c_U,...
    x_min, x_max, f_opt, x_opt);
Prob.P = P;
Prob.user = user;
% Setup function for the Trim Loss problems

function [A,b_L,b_U,c_L,c_U] = setup_trim(N,P,Nkmax,Bmax,delta,n,b,M,BIG)

if nargin < 9
    BIG = 1e6;
end

% Common matrices, zero and unit, resp.
ZP = spalloc(P,P,0);
IP = speye(P,P);

% (Bmax-delta)yj - sum_i(bi*rij) <=0, j=1,P
A1start = [ZP (Bmax-delta)*IP];
nrow = size((-b(1)*IP),1);
ncol = size((-b(1)*IP),2);
A1 = zeros(nrow,ncol*N);
for k=1:N
    A1(1:nrow,ncol*(k-1)+1:ncol*k) = -b(k)*IP;
end
A1 = [A1start A1];

% sum_i(bi*rij)-Bmax*yj <= 0, j=1,P
A2 = [ ZP -Bmax*IP -1*A1(:,2*P+1:end)];

% yj - sum_i(rij) <= 0, j=1,P
A3 = [ ZP IP          repmat(-IP,1,N) ];

% sum(rij)-Nkmax*yj <= 0, j=1,P
A4 = [ ZP   -Nkmax*IP repmat(IP,1,N)];

% yj-mj <= 0, j=1,P
A5 = [-IP         IP  repmat(ZP,1,N)];

% mj-Mj*yj <= 0, j=1,P
A6 = [IP -spdiags(M,0,P,P) repmat(ZP,1,N)];

% sum_j(mj) >= ...
A7 = [ones(1,P) zeros( 1, N*(P+1) ) ];

% y(k+1)-y(k) <= 0, k=1,2,...,P-1, NOTE P-1 rows!
nn = ones(P-1,1);
A8 = [ spalloc(P-1,P,0) spdiags([-nn nn],[0 1],P-1,P) spalloc(P-1,N*P,0) ];

% m(k+1)-m(k) <= 0, k=1,2,...,P-1, NOTE P-1 rows!
A9 = [ spdiags([-nn nn],[0 1],P-1,P) spalloc(P-1,(N+1)*P,0) ];

A = [A1;A2;A3;A4;A5;A6;A7;A8;A9];

c_L = n(:);
c_U = BIG*ones(N,1);

b_U = [     zeros(6*P,1) ; BIG ; zeros(2*(P-1),1) ];
b_L = [ -BIG*ones(6*P,1) ; ...
    max( ceil(sum(n)/Nkmax), ceil(sum(b.*n)/Bmax) ) ; ...
    -BIG*ones(2*(P-1),1) ];

% MODIFICATION LOG
%
% 021209 ango Wrote file
% 021227 hkh  Adding 9, 10, same as 3,1; Linear constraints as nonlinear
% 021229 hkh  Heavy modification of the problens
% 021230 hkh  Correcting with number in n(3), #8. Adding f_opt for #5-8
% 030124 ango Add x_min,x_max to all problems.
% 030206 ango Fix error in ConsPattern, problem 1
% 030206 hkh  Added comments to problem 1,2,9,10, Name change problem 9,10
% 030208 hkh  Added problem 11,12
% 041020 med  Added 13, 14
% 041117 med  xxx_prob removed and code added
% 050406 hkh  Use Inf instead of 1E6 for trim loss problems, avoid constraints
% 060927 med  Problem 14 x_0 set
% 070222 hkh  Change IntVars to use full vector format
% 080603 med  Switched to minlpAssign, cleaned
% 080916 nhq  Switched all x_opt into row-vectors
% 091008 hkh  Updated trim loss with correct optimal solutions from minlpSolve
