% Quality function from "A method for simulation based optimization using
% radial basis functions"
%
% On entry:
%    y           The point to evaluate the quality function in.
%    Prob        A problem structre where the following fields must be set:
%     .x_L       Lower bounds on variables.
%     .x_U       Upper bounds on variables.
%     .CGO
%      .X        The sampled points.
%      .min_sn   The minumum of the surrogate-model
%      .max_sn   Tha maximum of the surrogate-model
%      .QUAL
%       .sigma   The sigma parameter that balances local and global search.
%    
% sn_f is called from qual and its subroutines, so if using CGOLIB, the
% following also must be set:
%    Prob
%     .CGOLIB
%      .usecgolib  Must be set to 1 if using CGOLIB.
%      .rbf        The identifier of the RBF surface.
%
% This is just a model implementation. Much of this could and should be
% implemented in CGOLIB. Issues:
%
%  * How to maximize the quality function?
%    As we are using Monte Carlo Integration, we need a whole lot of calls
%    to the integrand, which itself contains calls to sn_f and computations
%    of distances. We have no analytical derivative, and because of the huge
%    amounts of calls that are needed in order to get good accuracy in the
%    integral, it is probably not suitable for numerical differentiation
%    either. Well, so maybe DIRECT, or even another black-box CGO solver
%    would be sutiable?
%
%    (Hum, thinking. Maybe we could find an analytical derivative of the
%    integrand and integrate it in order to get a derivative? A problem is 
%    that the non-zero region of the integral domain is moving as y changes.) 
%
%  * How to minimize the computation time?
%    1. Find an efficient way of computing the "smallest enclosing box" of the
%       non-zero simplex. The two optimzation problems per dimension are
%       costly.
%    2. Find a way of removing redundant constraints. In this process, we
%       would automatically get information of what points are candidates
%       of being the closest point to any point within the simplex. This
%       would decrease the time needed to compute the norm in qualInt,
%       because we would have to loop over a less number of points.
%    3. Is there any other way to integrate the function than Monte
%       Carlo Integration?
%   
% The text above was written 080729 by Fredrik Hellman
%
% After some testing, several interesting things have turned out:
%
%  1. The enclosing box doesn't make sense for higher dimensions.
%     Even though we have found the smalles enclosing box, the
%     ratio: simplex volume over box volume becomes smaller and smaller
%     as the dimension increase. So the enclosing box doesn't work
%     for higher dimensions. We need something else to make the
%     integration. Maybe: Generate Gaussian-points within the simplex
%     somehow. One method is here: 
%     
%     http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=9435&objectType=File
%
%     But it requires the vertices of the simplex, and we only have
%     some inequalities. In that case, we would like to find all basic
%     feasible solutions to an LP-problem with our inequalities as
%     constraints. Can this be done efficiently? Found this from 1961:
%
%      An Algorithm for Finding All Vertices of Convex Polyhedral Sets
%      M. L. Balinski
%
%  2. In the code is now implemented a way of removing redundant
%     constraints. It is based on computing the convex hull of a polyhedron
%     (the dual to the simplex). This doesn't work well for higher
%     dimensions and doesn't cause any speed up.
%
%  3. A method for determining the absolute tolerance of the quadrature was
%     implemented. However, it was based on computing the volume of the
%     simplex to get an estimate of the integrand. Computing the volume was
%     a convex hull-task and hence really slow. It is out-iffed (if 0)
%     below.
%
% The text above was written 080801 by Fredrik Hellman

function [q] = qual(y, Prob)

% Set some options. They are hardcoded now, but could be set through
% some option structure in the future.

% BOXMAXDIM - Don't use enclosing box for dimension greater than
%             BOXMAXDIM.
BOXMAXDIM = 20;

% NREDUNDMAXDIM - Don't remove redundant constraints for dimension
%                 greater than NOREDUNDMAXDIM.

NOREDUNDMAXDIM = 5;

% Further below, we also set absolute tolerance for quadrature and number
% of sample points for the Monte Carlo integration. This must be
% automatically adjusted somehow.

d = length(y);

X = Prob.CGO.X;
n = size(X,2);

% Find the simplex
[A, b_U] = clsmplx(X, y);

% Add box bounds
A = [A; -eye(d); eye(d)];
b_U = [b_U; Prob.x_L; Prob.x_U];

% Remove redundant constraints
if(d > 1 && d <= NOREDUNDMAXDIM)
    [An, b_Un, nr] = noredund(A, b_U, y);
else
    An = A;
    b_Un = b_U;
    nr = 1:n;
end

% Some of the box bound constraints may be are removed.
nr(nr>n) = [];

Prob.CGO.QUAL.y = y;
Prob.CGO.QUAL.A = An;
Prob.CGO.QUAL.b_U = b_Un;
Prob.CGO.QUAL.nr = nr;

% Find the enclosing box
% For bigger dimensions, using the enclosing box has no
% significant effect on the performance, because the box
% will basically be the original box bounds.

x_L = Prob.x_L;
x_U = Prob.x_U;

if(d <= BOXMAXDIM)
    SolverLP = GetSolver('lp');

    c = zeros(d,1);
    BoxProb = lpAssign(c, An, [], b_Un, Prob.x_L, Prob.x_U, y);

    for i=1:d
        BoxProb.QP.c = c;
        BoxProb.QP.c(i) = 1;
        R = tomRun(SolverLP,BoxProb);
        x_L(i) = R.x_k(i);

        BoxProb.QP.c(i) = -1;
        R = tomRun(SolverLP,BoxProb);
        x_U(i) = R.x_k(i);
    end
end

fprintf('Box:\n');
disp([x_L, x_U]);

% To determine the tolerance of the integral calculation,
% an estimate of the integral is computed.
% The largest possible value of the integrand is found,
% and multiplied by the volume of the simplex to
% obtain an upper limit of the value of the integral.
% We can take this as our estimate.
% Then we multiply this by the relative tolerance we
% want.
%
% The volume of the simplex is somewhat difficult to
% get. In noredund, we have the vertices of the dual
% polyhedron. I would like to use these to get
% the volume of the polyhedron, but I haven't found
% an easy way. Instead I found some code:
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=7894&objectType=FILE
% which give me the vertices of our simplex.
% Then I use convhulln to get the volume. This is 
% probably a really inefficient way of doing this...

% This doesn't work very well. The estimate must
% depend on the dimension as well somehow...
% This method doesn't only give a bad estimate,
% it is also inefficient.
if 0
    Anv = con2vert(An, b_Un, y);
    [crap, V] = convhulln(Anv);
    
    if d==2
        figure
        hold on;
        plot(Anv(:,1), Anv(:,2),'.');
        plot(Prob.CGO.X(1,:), Prob.CGO.X(2,:), 'r.');
        plot(y(1,:), y(2,:), 'ks');
        pause
    end
    
    % The integrand is at most equal to
    %  U_X - U_Xy
    % which is at most min_(w in X) ||y-w||.
    
    maxInt = inf;
    for i=1:length(nr)
        dY = norm(y-Prob.CGO.X(:,nr(i)));
        if(maxInt > dY)
            maxInt = dY;
        end
    end
    
    % omega is at most 1.
    
    % maxq = V*maxInt; % NOT USED
end

% Integrate
% =========

% We need a way to set the tolerance for quadrature!
intTol = 1e-7;

% We also need a way to set the number of monte carlo
% points.
N = 10000;

monte = 0;
if(d == 1 && ~monte)
    q = quad(@qualInt, x_L, x_U, intTol, 0, Prob);
elseif(d == 2 && ~monte)
    q = dblquad(@qualIntDbl, x_L(1), x_U(1), x_L(2), x_U(2), intTol, [], Prob);
elseif(d == 3 && ~monte)
    q = triplequad(@qualIntTriple, x_L(1), x_U(1), x_L(2), x_U(2), x_L(3), x_U(3), intTol, [], Prob);
else
    qs = 0;
    % bV = Box volume
    bV = prod(x_U-x_L);

    % nz - number of nonzero integrands.
    nz = 0;
    for i=1:N
        x = x_L + rand(d,1).*(x_U-x_L);
        qi = qualInt(x, Prob);
        qs = qs + qi;
        if(qi > 0)
            nz = nz + 1;
        end
    end
    q = bV*qs/N;

    % nz/N gives us the hit ratio of the samples.
    fprintf('MC hit ratio: ');
    disp(nz/N);
    
    % 1000 is no magic number. I just picked something to enlight the problem!
    if nz < 1000
        warning('Less than 1000 points were sampled within the simplex');
    end
end


% qualInt
% The integrand.
%
% Prob.CGO.QUAL.
%  A
%  b_U     A and b_U form the simplex where the integrand is nonzero
%  nr      These are the indices of the points forming the simplex.
%  y
%
% Prob.CGO.
%  X
%  min_sn  Minimum of s_n
%  max_sn  Maximum of s_n

function [f] = qualInt(xv, Prob)

f = zeros(size(xv,2),1);

for j = 1:size(xv,2)
    x = xv(:,j);
    
    if(all(Prob.CGO.QUAL.A*x <= Prob.CGO.QUAL.b_U))
    
        y = Prob.CGO.QUAL.y;
        min_sn = Prob.CGO.min_sn;
        max_sn = Prob.CGO.max_sn;
        nr = Prob.CGO.QUAL.nr;
    
        % Compute minimum distance between x and the points in X (only the
        % nr-points, those forming the simplex)
        U_X = inf;
        
        if 1
            for i=1:length(nr)
                dX = norm(x-Prob.CGO.X(:,nr(i)));
                if(U_X > dX)
                    U_X = dX;
                end
            end
        else % Here is the old version, comparing distance to all
             % points.
            for i=1:size(Prob.CGO.X,2)
                dX = norm(x-Prob.CGO.X(:,i));
                if(U_X > dX)
                    U_X = dX;
                end
            end
        end
        
        % Minimum distance between x and (X union {y}) is the distance between x and y
        U_Xy = norm(x-y);
    
        Prob.CGO.fnStar = NaN;
        Prob.CGOLIB.usecgolib = 1;
        sn_x = sn_f(x, Prob);
    
        omega = exp(-Prob.CGO.QUAL.sigma*(sn_x-min_sn)/(max_sn-min_sn));
        omega = min(1, max(omega, 0));
    
        f(j) = (U_X - U_Xy)*omega;
    else
        f(j) = 0;
    end
end

function f = qualIntDbl(xv1, x2, Prob)
    xv = [xv1; repmat(x2, 1, size(xv1,2))];
    f = qualInt(xv, Prob);
    

function f = qualIntTriple(xv1, x2, x3, Prob)
    xv = [xv1; repmat(x2, 1, size(xv1,2)); repmat(x3, 1, size(xv1,2))];
    f = qualInt(xv, Prob);
    
% clsmplx
% Find the simplex in which all points are closer to y than to any other point
% in X. Returns A and b_U in A*x <= b_U where x are the points in the
% simplex. Nothing is done to remove redundant inequalities.

function [A, b_U] = clsmplx(X, y)

n = size(X,2);
d = size(X,1);
A = zeros(n,d);
b_U = zeros(n,1);

for i=1:n
    d = X(:,i)-y(:);
    b = d'*X(:,i)-d'*d/2;
    A(i,:) = d';
    b_U(i) = b;
end
