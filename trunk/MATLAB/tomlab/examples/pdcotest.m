function pdcotest( m,n )

% m=50;  n=100;  pdcotest( m,n );
% Defines an m by n LP problem
% and runs it on pdco.m with the entropy objective function.

if nargin < 1
   help pdcotest
   return
end

%-----------------------------------------------------------------------
% 08 Oct 2002: Simple test program for pdco.m.
%              "A" is an explicit sparse matrix, not a function.
%              Michael Saunders, SOL, Stanford University.
%-----------------------------------------------------------------------

  [A,b,bl,bu,d1,d2] = toydata( m,n );   % Private function below
  D  = sum(A,1);   D(find(D==0)) = 1;
  D  = sparse( 1:n, 1:n, 1./D, n, n );
  A  = A*D;                             % Normalize cols of A

  options = pdcoSet;

  x0 = ones(n,1)/n;   xsize = 1/n;
  y0 = zeros(m,1);
  z0 = ones(n,1)/n;   zsize = 1/n;

  options.mu0       = 1e-1;  % 1.0 starts near central path
  options.LSQRatol1 = 1e-6;  % Let LSQR solve loosely to start with
  options.wait      = 1;     % Allow options to be reviewed before solve
  options.Method    = 4;     % Use Tomlab Special PDCO LSQR MEX file 
  Prob.PriLevOpt    = 1;

  [x,y,z,inform,PDitns,CGitns,time] = ...
    pdco( 'entropy',A,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize,Prob );

  %plot(b-A*x)

% keyboard                   % Allow review of x,y,z, etc.
%-----------------------------------------------------------------------
% End function pdcotest
%-----------------------------------------------------------------------


function [A,b,bl,bu,d1,d2] = toydata( m,n )

%        [A,b,bl,bu,d1,d2] = toydata( m,n );
%        defines an m by n matrix A and rhs vector b,
%        for use with pdco.m.

%-----------------------------------------------------------------------
% 12 Feb 2001: First version of toydata.m.
% 30 Sep 2002: pdsco version modified for pdco.
%-----------------------------------------------------------------------

  rand('state',0);
  density = 0.50;
  rc      = 1e-1;

  em      = ones(m,1);
  en      = ones(n,1);
  zn      = zeros(n,1);
  bigbnd  = 1e+30;

  if m==1
    A = sparse(ones(m,n));
  else
    A = [sprand(m-1,n,density,rc)
         en'                     ];
  end
  x       = en/n;
  b       = full(A*x);

  bl      = zn;
  bu      = en;
% bl(1)   = - bigbnd;   % Test "free variable" with no bounds
% bu(1)   =   bigbnd;
  gamma   = 1e-4;
  delta   = 1;          % Least squares
  d1      = en*gamma;
  d2      = em*delta;
  d2(m)   = 1e-4;       % Make e'x = 1 satisfied more accurately.

%-----------------------------------------------------------------------
% End private function toydata
%-----------------------------------------------------------------------
