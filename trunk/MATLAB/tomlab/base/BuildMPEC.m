% function MProb = BuildMPEC(Prob,mpec)
%
% Converts a regular optimization problem in the TOMLAB format to a
% complementarity (equilibrium) problem by adding slack variables and
% modifying the constraint structure if necessary.
%
% MPEC (Mathematical Program with Equlibrium Constraints) contain pairs
% of variables and/or constraints which are designated to be complementary, e.g.
%
%  0 <= x(i) _|_ x(2) >= 0
%
% a condition that implies that either (or both) of x(1),x(2) must be zero.
%
% This can be extended to constraints as well - suppose a problem has two
% nonlinear constraints c1(x) and c2(x):
%
% 0 <= c1(x)
% 0 <= c2(x)
%
% By adding two slack variables to the problem and specifying these as
% complementary, the constraints will in fact become complementary:
%
% c1(x) - s(1) == 0
% c2(x) - s(2) == 0
% 0 <= s(1) _|_ s(2) >= 0
%
% The complementary pairs of variables and constraints are described as rows
% in the input argument 'mpec', a matrix with 6 columns and any number of rows.
%
% Each row must have exactly two nonzero elements describing which
% variable(s) and/or constraint(s) are part of the complementarity pair.
%
% Columns 1:2 refer to variables, columns 3:4 to linear constraints,
% and 5:6 to nonlinear constraints:
%
% >> mpec = [   var1,var2 , lin1,lin2 , non1,non2  ; ... ];
%
% To specify more than one pair, simply add more rows.
%
% EXAMPLE:
%
% To make variables 1 and 2 complementary 0 <= x(1) _|_ x(2) >= 0, do:
%
% >> mpec = [ 1,2, 0,0, 0,0 ];
%
% >> Prob2 = BuildMPEC(Prob,mpec);
%
% If only existing variables are chosen for complementarity pairs, no slack
% variables will be added. If constraints are included in a complementarity
% pair, one slack variable will be added for each constraint.
%
% EXAMPLE:  Suppose the first linear and the first nonlinear constraint
% should be complementary:    b_L <= A(1,:)*x _|_ c1(x) >= 0
%
% The mpec row can be specified as [ 0,0 , 1,0 , 1,0 ].
%
% BuildMPEC will add slacks and insert gateway routines into Prob.FUNCS that
% automatically handles the new slack variables. The linear
% constraint matrix Prob.A will also be expanded with the slack
% coefficient.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written May 19, 2006.   Last modified Aug 14, 2006.

function MProb = BuildMPEC(Prob,mpec)

if nargin < 2
    error('BuildMPEC needs two input arguments, the Prob structure and the list of complementarity pairs mpec')
end

if isempty(mpec)
    if(Prob.Warning)
        warning('mpec input argument in BuildMPEC is empty, returning Prob without changes');
        MProb = Prob;
        return;
    end
end

MProb = Prob;

% x_L = Prob.x_L;
% x_U = Prob.x_U;

A   = Prob.A;
b_L = Prob.b_L;
b_U = Prob.b_U;

c_L = Prob.c_L;
c_U = Prob.c_U;

ConsPattern = DefPar(Prob,'ConsPattern',[]);

BIG = DefPar(Prob,'BIG',1e20);

% Number of requested mpec pairs
nmpec = size(mpec,1);

% Current problem dimensions
N  = Prob.N;
m1 = Prob.mLin;
m2 = Prob.mNonLin;

% MP is the extension to dc/ConsPattern - if any constraints are
% to be included in the complementarity

MP  = sparse(m2,nnz(mpec(:,3:6)));
mp  = [];
ns  = 0;
exA = [];
exC = [];

for k=1:nmpec

    if nnz(mpec(k,:)) > 2
        error(['Too many nonzeros in mpec ' num2str(k) ]);
    elseif nnz(mpec(k,:)) < 2
        error(['Too few nonzeros in mpec ' num2str(k) ]);
    end

    x1i = mpec(k,1);
    x2i = mpec(k,2);

    a1i = mpec(k,3);
    a2i = mpec(k,4);

    c1i = mpec(k,5);
    c2i = mpec(k,6);

    i1 = max(x1i,x2i);
    i2 = max(a1i,a2i);
    i3 = max(c1i,c2i);

    % Check validity of pair indices
    if(i1>N)
        e = sprintf('mpec pair %d includes variable %d which is higher than Prob.N = %d',k,i1,N);
        error(e);
    end

    if(i2>m1)
        e = sprintf('mpec pair %d includes linear constraint %d which is higher than Prob.mLin = %d',k,i2,m1);
        error(e);
    end

    if(i3>m2)
        e = sprintf('mpec pair %d include nonlinear constraint %d which is higher than Prob.mNonLin = %d',k,i3,m2);
        error(e);
    end

    if( all([c1i,c2i]) )    % Two nonlinear constraints - add two slacks
        [ConsPattern,c_L,c_U,k1,exC] = XNonLinear(ConsPattern,c_L,c_U,c1i,0,exC,BIG);
        [ConsPattern,c_L,c_U,k2,exC] = XNonLinear(ConsPattern,c_L,c_U,c2i,0,exC,BIG);

        MP(k1,ns+1) = 1;
        MP(k2,ns+2) = 1;

        % Enter new slacks as a mpec pair
        mp = [mp ; N+ns+1, N+ns+2 ];
        ns = ns + 2;

    elseif ( all([a1i,a2i]) ) % Two linear constraints - add two slacks

        %       A(a1i,N+ns+1) = -1.0;
        %       A(a2i,N+ns+2) = -1.0;
        %       b_U(a1i) = b_L(a1i);
        %       b_U(a2i) = b_L(a2i);
        %       mp = [mp ; N+ns+1, N+ns+2];
        %       ns = ns + 2;

        [ A, b_L, b_U, k1, exA ] = XLinear( A, b_L, b_U, a1i, 0, exA, BIG );
        [ A, b_L, b_U, k2, exA ] = XLinear( A, b_L, b_U, a2i, 0, exA, BIG );

        A(k1,N+ns+1) = -1.0;
        A(k2,N+ns+2) = -1.0;

        % Enter new slacks as a mpec pair
        mp = [mp ; N+ns+1, N+ns+2];
        ns = ns + 2;

    elseif( all([x1i,x2i]) ) % Two variables - no slacks
        mp = [mp ; x1i,x2i];

    elseif( i1~=0 & any([a1i,a2i]) ) % One linear constraint and one existing variable (i1)
        [ A, b_L, b_U, k, exA ] = XLinear( A, b_L, b_U, a1i, a2i, exA, BIG );

        A(k,N+ns+1) = -1.0;
        mp = [mp ; i1,N+ns+1];
        ns = ns + 1;

    elseif( i1>0 & any([c1i,c2i]) ) % One nonlinear constraint and one existing variable.
        [ConsPattern,c_L,c_U,k,exC] = XNonLinear(ConsPattern,c_L,c_U,c1i,c2i,exC,BIG);

        MP( k, ns+1 ) = 1;
        mp = [mp ; N+ns+1,i1];
        ns = ns + 1;

    elseif( any([a1i,a2i]) & any([c1i,c2i]) ) % One linear and one nonlinear constraint
        [ConsPattern,c_L,c_U,k,exC] = XNonLinear(ConsPattern,c_L,c_U,c1i,c2i,exC,BIG);
        MP( k, ns+1 ) = 1;

        [A, b_L, b_U, k, exA ] = XLinear( A, b_L, b_U, a1i, a2i, exA, BIG);
        A( k, N+ns+2 ) = -1.0;

        mp = [mp ; N+ns+1,N+ns+2];
        ns = ns + 2;
    end
end

MPEC = [];

MPEC.mp = mp;
MPEC.MP = MP;
MPEC.ns = ns;
MPEC.exC = exC;
MPEC.exA = exA;

MProb.MPEC    = MPEC;
MProb.orgProb = Prob;

MProb.KNITRO.mpec = mp;

% New problem has N+sn variables
MProb.N = N + ns;

% Slacks are all 0.0 <= s <= inf
MProb.x_L(N+1:N+ns) = 0.0;
MProb.x_U(N+1:N+ns) = Inf;

MProb.b_L = b_L;
MProb.b_U = b_U;

MProb.c_L = c_L;
MProb.c_U = c_U;

% Starting values for slacks is zero (set if a starting point exists)
if ~isempty(MProb.x_0)
    MProb.x_0(N+1:N+ns)=0.0;
end

%
% Replace user functions if slacks have been added.
if ns>0
    if ~isempty(MProb.FUNCS.c)
        MProb.FUNCS.c = 'mpec_c';
    end

    if ~isempty(MProb.FUNCS.dc)
        MProb.FUNCS.dc = 'mpec_dc';
    elseif MProb.CheckNaN>0
         % No dc(x) is given, but CheckNaN request. mpec_dc can use this.
         MProb.FUNCS.dc = 'mpec_dc';
    end

    if ~isempty(MProb.FUNCS.d2c)
        MProb.FUNCS.d2c = 'mpec_d2c';
    end

    if ~isempty(MProb.FUNCS.f)
        MProb.FUNCS.f = 'mpec_f';
    end

    if ~isempty(MProb.FUNCS.g)
        MProb.FUNCS.g = 'mpec_g';
    end

    if ~isempty(MProb.FUNCS.H)
        MProb.FUNCS.H = 'mpec_H';
    end
end

%
% Expand ConsPattern with MP nonzero structure
if ~isempty(ConsPattern) & ns>0
    if isempty(MP)
        ConsPattern(end,N+ns) = 0; %  = [ ConsPattern sparse(m2,ns) ];
    else
        ConsPattern(1:size(MP,1),N+1:N+ns) = MP; %  = [ ConsPattern spones(MP) ];
    end
end
MProb.ConsPattern = ConsPattern;
MProb.mNonLin = length(c_L);

% Expand relevant patterns with rows and/or columns for slacks. All slacks
% enter linearly so only zero columns need to be added.

if m1>0 & ns>0 & size(A,2)<N+ns
    A(end,N+ns) = 0.0;
end

MProb.A    = A;
MProb.mLin = size(A,1);

if ~isempty(MProb.HessPattern) & ns > 0
    MProb.HessPattern(N+ns,N+ns) = 0.0;
end

if ~isempty(MProb.d2LPattern) & ns > 0
    MProb.d2LPattern(N+ns,N+ns) = 0.0;
end

if ~isempty(MProb.d2cPattern) & ns > 0
    MProb.d2cPattern(N+ns,N+ns) = 0.0;
end

return
% End of BuildMPEC main function


% XNonLinear

function [ConsPattern,c_L,c_U,k,exC] = XNonLinear(ConsPattern,c_L,c_U,c1i,c2i,exC,BIG)

k = max(abs( [c1i,c2i] ));

% Analyze bounds. Two more or less "wrong" cases exist:
% - if equality, change nothing.
% - if free, throw an error because KNITRO will fail on this anyway.

if c_U(k) == c_L(k)
    s = sprintf('Nonlinear constraint %d is already an equality',k);
    warning(s);
elseif c_L(k) <= -BIG & c_U(k) >= BIG % Free row - KNITRO will not cope, thus an error.
    s = sprintf('Nonlinear constraint %d is free',k);
    error(s);
elseif c_L(k) <= -BIG & c_U(k) < BIG % Only upper bounded
    % Use upper bound in MPEC
    c_L(k) = c_U(k);
elseif c_L(k) > -BIG & c_U(k) >= BIG % Only lower bounded
    % Use lower bound in MPEC
    c_U(k) = c_L(k);
else
    % Bounded both upper and lower - use sign of mpec pair indicator to
    % determine which bound to use. Lower if negative.
    % Also, duplicate the "other half" of the constraint.

    if min(c1i,c2i) < 0 % Lower half becomes MPEC constraint
        % Duplicate upper half of original constraint

        if ~isempty(ConsPattern)
            ConsPattern(end+1,:) = ConsPattern(k,:);
        end

        c_L(end+1) = -Inf;
        c_U(end+1) = c_U(k);

        c_U(k) = c_L(k);
    else
        % Duplicate lower half
        if ~isempty(ConsPattern)
            ConsPattern(end+1,:) = ConsPattern(k,:);
        end

        c_L(end+1) = c_L(k);
        c_U(end+1) = Inf;

        c_L(k)     = c_U(k);
    end

    % Keep track of where the extra constraint came from
    exC = [ exC ; k];
end

% end of XNonLinear

% XLinear

function [A,b_L,b_U,k,exA] = XLinear(A,b_L,b_U,a1i,a2i,exA,BIG)

k = max(abs( [a1i,a2i] ));

% Analyze bounds. Two more or less "wrong" cases exist:
% - if equality, change nothing.
% - if free, throw an error because KNITRO will fail on this anyway.

if b_U(k) == b_L(k)
    s = sprintf('Linear constraint %d is already an equality constraint',k);
    warning(s);
elseif b_L(k) <= -BIG & b_U(k) >= BIG % Free row - KNITRO will not cope, thus an error.
    s = sprintf('Linear constraint %d is free',k);
    error(s);
elseif b_L(k) <= -BIG & b_U(k) < BIG % Only upper bounded
    % Use upper bound in MPEC
    b_L(k) = b_U(k);
elseif b_L(k) > -BIG & b_U(k) >= BIG % Only lower bounded
    % Use lower bound in MPEC
    b_U(k) = b_L(k);
else
    % Bounded both upper and lower - use sign of mpec pair indicator to
    % determine which bound to use. Lower if negative.
    % Also, duplicate the "other half" of the constraint.

    if min(a1i,a2i) < 0 % Lower part of A(k,:) becomes MPEC constraint
        % Duplicate upper half of original constraint
        A(end+1,:) = A(k,:);
        b_L(end+1) = -Inf;
        b_U(end+1) = b_U(k);

        b_U(k) = b_L(k);
    else
        % Duplicate lower half
        A(end+1,:) = A(k,:);
        b_L(end+1) = b_L(k);
        b_U(end+1) = Inf;
        b_L(k)     = b_U(k);
    end

    % Keep track of where the extra constraint came from
    exA = [ exA ; k];
end

% end of XLinear

% MODIFICATION LOG
%
% 060519 ango Wrote file
% 060524 med  Minor updates, x_L, x_U and W commented
% 060814 med  FUNCS used for callbacks instead
