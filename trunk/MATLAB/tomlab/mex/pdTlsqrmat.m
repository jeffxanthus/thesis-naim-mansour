% pdTlsqrmat is required by pdco.m (when it calls pdTlsqr.m).
% It forms Mx or M'x for some operator M that depends on LSproblem below.
%
% mlsqr, nlsqr  are the dimensions of the LS problem that lsqr is solving.
%
% Aname is pdsco's Aname.
%
% rw contains parameters [explicit LSproblem LSmethod LSdamp]
% from pdco.m to say which least-squares subproblem is being solved.
%
% global pdDDD1 pdDDD3 provides various diagonal matrices
% for each value of LSproblem.

%-----------------------------------------------------------------------
% 17 Mar 1998: First version to go with pdsco.m and lsqr.m.
% 01 Apr 1998: global pdDDD1 pdDDD3 now used for efficiency.
% 11 Feb 2000: Added diagonal preconditioning for LSQR, LSproblem = 1.
% 14 Dec 2000: Added diagonal preconditioning for LSQR, LSproblem = 12.
% 30 Jan 2001: Added diagonal preconditioning for LSQR, LSproblem = 21.
% 12 Feb 2001: Included in pdco.m as private function.
%              Specialized to allow only LSproblem = 1.
% 03 Oct 2002: First version to go with pdco.m with general H2 and D2.
% 16 Oct 2002: Aname is now the user's Aname.
%
% 20 Jan 2003: Adapted to Tomlab by Kenneth Holmstrom
%-----------------------------------------------------------------------

function y = pdTlsqrmat( mode, mlsqr, nlsqr, x, Aname, rw )

global pdDDD1 pdDDD2 pdDDD3

LSproblem = rw(2);
precon    = logical(rw(7));

if LSproblem == 1
    % The operator M is [ D1 A'; D2 ].
    m = nlsqr;
    n = mlsqr - m;
    if mode == 1
        if precon, x = pdDDD3.*x; end
        t = pdTmat( Aname, 2, m, n, x );   % Ask 'aprod' to form t = A'x.
        %t = Aname'*x;
        y = [ (pdDDD1.*t); (pdDDD2.*x) ];
    else
        t = pdDDD1.*x(1:n);
        y = pdTmat( Aname, 1, m, n, t );   % Ask 'aprod' to form y = A t.
        %y = Aname*t;
        y = y + pdDDD2.*x(n+1:mlsqr);
        if precon, y = pdDDD3.*y; end
    end
else
    error('Error in pdTlsqrmat: Only LSproblem = 1 is allowed')
end

function y = pdTmat( Aname, mode, m, n, x )

%        y = pdTmat( Aname, mode, m, n, x )
%    computes y = Ax (mode=1) or A'x (mode=2)
%    for a matrix A defined by pdco's input parameter Aname.

%-----------------------------------------------------------------------
% 04 Apr 1998: Default A*x and A'*y function for pdco.m.
%              Assumed A was a global matrix pdAAA created by pdco.m
%              from the user's input parameter A.
% 16 Oct 2002: pdAAA eliminated to save storage.
%              User's parameter Aname is now passed thru to here.
% 01 Nov 2002: Bug: feval had one too many parameters.
%
% 20 Jan 2003: Adapted to Tomlab by Kenneth Holmstrom
%-----------------------------------------------------------------------

if ischar(Aname)
    y = feval( Aname, mode, m, n, x );
else
    if mode==1,  y = Aname*x;  else  y = Aname'*x;  end
end

%-----------------------------------------------------------------------
% End private function pdTmat
%-----------------------------------------------------------------------