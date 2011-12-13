function pdcotest4( m,n )

% m=50;  n=100;  pdcotest4( m,n );
% Defines an m by n LP problem
% and runs it on pdco.m with local linear objective function,
% first with loose opttol, then again to find a correction.

if nargin < 1
   help pdcotest4
   return
end

%-----------------------------------------------------------------------
% 23 Oct 2002: Simple test program for pdco.m.
%              "A" is an explicit sparse matrix, not a function.
%              "phi(x)" is defined by the private function myOBJ(x);
%              Michael Saunders, SOL, Stanford University.
%-----------------------------------------------------------------------

  global pdOBJVEC pdOBJADD
  
  [A,b,bl,bu,c,d1,d2] = toydata( m,n );   % Private function below
  D  = sum(A,1);   D(find(D==0)) = 1;
  D  = sparse( 1:n, 1:n, 1./D, n, n );
  A  = A*D;                               % Normalize cols of A

  pdOBJ    = @myOBJ;      % Test function_handle
  pdOBJVEC = c;           % Inelegant way to make c available to myOBJ
  pdOBJADD = 0;
  
  options = pdcoSet;

  en = ones(n,1);
  zn = zeros(n,1);
  zm = zeros(m,1);
  x0 = en/n;   xsize = 1;
  y0 = zm;
  z0 = en/n;   zsize = 1;

  options.mu0       = 1e-1;  % 1.0 starts near central path
  options.LSQRatol1 = 1e-6;  % Let LSQR solve loosely to start with
  options.wait      = 1;     % Allow options to be reviewed before solve

  options.FeaTol    = 1e-3;  % Loose solve
  options.OptTol    = 1e-3;  % Loose solve
  Prob.PriLevOpt    = 1;
  
  [x1,y1,z1,inform,PDitns,CGitns,time] = ...
    pdco( pdOBJ,A,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize, Prob );

  bl1 = bl - x1;
  bu1 = bu - x1;
  b1  = b  - A*x1;
  x0  = zn;  y0 = zm;  z0 = zn;  % New solution is approx x=0, y=0, z=0

  d   = A'*y1 + z1;
  pdOBJVEC = c - d;              % Make new c available to myOBJ
  pdOBJADD = c'*x1;

  options.x0min     = 1e-2;
  options.z0min     = 1e-2;
  options.LSQRatol1 = options.LSQRatol1*1e-2;
  options.Method    = 4;         % Use Tomlab Special PDCO LSQR MEX file 

% xsize = options.OptTol * 10;   % Fairly small
  zsize = norm(z1,inf);          % z1 should be about right
  xsize = 1;
  zsize = 1;

  options.OptTol    = 1e-3;      % HOPE LOOSE TOL IS GOOD ENOUGH!
  options.OptTol    = 1e-6;      % Try tight one anyway
  
  [x2,y2,z2,inform,PDitns,CGitns,time] = ...
    pdco( pdOBJ,A,b1,bl1,bu1,d1,d2,options,x0,y0,z0,xsize,zsize,Prob );

  x     = x1+x2;
  y     = y1+y2;
  z     = z1+z2;
  rnorm = norm(A*x-b,inf)
  dnorm = norm(A'*y+z-c,inf)
  %keyboard                   % Allow review of x,y,z, etc.
%-----------------------------------------------------------------------
% End function pdcotest4
%-----------------------------------------------------------------------


function [A,b,bl,bu,c,d1,d2] = toydata( m,n )

% Defines an m by n matrix A and associated vectors
% for testing pdco.m on a linear program.

%-----------------------------------------------------------------------
% 23 Oct 2002: Entropy function in pdcotest.m makes problem "too easy".
%              Linear programs are harder for primal-dual barrier.
%-----------------------------------------------------------------------

  rand('state',0);
  density = 0.50;
  rc      = 1e-1;

  em      = ones(m,1);
  en      = ones(n,1);
  zn      = zeros(n,1);

  A       = sprand(m,n,density,rc);
  x       = en;
  b       = full(A*x);
  c       = rand(n,1);

  bl      = - 10*en;
  bu      =   10*en;
% bl(1)   = -inf;       % Test "free variable" with no bounds
% bu(1)   = +inf;
  d1      = 1e-4;       % pdco will use D1 = d1*I
  d2      = 1e-3;       % pdco will use D2 = d2*I

%-----------------------------------------------------------------------
% End private function toydata
%-----------------------------------------------------------------------


function [obj,grad,hess] = myOBJ( x, Prob )
%        [obj,grad,hess] = myOBJ( x, Prob )
%        computes the objective value, gradient and diagonal Hessian
%        of a separable convex function, for use with pdco.m.
%        This is an example LINEAR objective function.
%        It tests use of a function_handle to find this private function.
%        Same as linobj.m but picks up c from global pdOBJVEC.
%        Also adds global constant pdOBJADD to obj.

  global pdOBJVEC pdOBJADD

  n    = length(x);
  c    = pdOBJVEC;
  obj  = c'*x + pdOBJADD;
  grad = c;
  hess = zeros(n,1);

%-----------------------------------------------------------------------
% End private function myOBJ
%-----------------------------------------------------------------------
