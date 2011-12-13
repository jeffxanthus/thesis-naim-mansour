% TLSQR finds a solution x to the following problems:
%
% 1. Unsymmetric equations - solve  A*x = b
%
% 2. Linear least squares  - solve  A*x = b in the least-squares sense
%
% 3. Damped least squares  - solve  (   A    )*x = ( b )
%                                   ( damp*I )     ( 0 ) in least-squares sense
%
% where A is a matrix with m rows and n columns, b is an m-vector, and damp is
% a scalar (All quantities are real). The matrix A is intended to be a large
% and sparse array, but both full and sparse Matlab arrays are handled.
%
% Note:  LSQR uses an iterative method to approximate the solution.
% The number of iterations required to reach a certain accuracy
% depends strongly on the scaling of the problem.  Poor scaling of
% the rows or columns of A should therefore be avoided where possible.
%
% For example, in problem 1 the solution is unaltered by
% row-scaling.  If a row of A is very small or large compared to the other
% rows of A, the corresponding row of ( A  b ) should be scaled up or down.
%
% In problems 1 and 2, the solution x is easily recovered following
% column-scaling.  Unless better information is known, the nonzero columns of
% A should be scaled so that they all have the same Euclidean norm (e.g., 1.0).
%
% In problem 3, there is no freedom to re-scale if damp is nonzero. However,
% the value of damp should be assigned only after attention has been paid to
% the scaling of A.
%
% The parameter damp is intended to help regularize
% ill-conditioned systems, by preventing the true solution from
% being very large.  Another aid to regularization is provided by
% the parameter aCond, which may be used to terminate iterations
% before the computed solution becomes very large.
%
% Note that x is not an input parameter.  If some initial estimate x0 is known
% and if damp = 0, one could proceed as follows:
%
%   1. Compute a residual vector     r0 = b - A*x0.
%   2. Use LSQR to solve the system  A*dx = r0.
%   3. Add the correction dx to obtain a final solution x = x0 + dx.
%
% This requires that x0 be available before and after the call
% to LSQR.  To judge the benefits, suppose LSQR takes k1 iterations
% to solve A*x = b and k2 iterations to solve A*dx = r0.
% If x0 is "good", norm(r0) will be smaller than norm(b).
% If the same stopping tolerances aTol and bTol are used for each
% system, k1 and k2 will be similar, but the final solution x0 + dx
% should be more accurate.  The only way to reduce the total work
% is to use a larger stopping tolerance for the second system.
% If some value bTol is suitable for A*x = b, the larger value
% bTol*norm(b)/norm(r0)  should be suitable for A*dx = r0.
%
% Preconditioning is another way to reduce the number of iterations.
% If it is possible to solve a related system M*x = b efficiently,
% where M approximates A in some helpful way
% (e.g. M - A has low rank or its elements are small relative to
% those of A), LSQR may converge more rapidly on the system
%       A*M(inverse)*z = b,
% after which x can be recovered by solving M*x = z.
%
% NOTE: If A is symmetric, LSQR should not be used!
% Alternatives are the symmetric conjugate-gradient method (cg) and/or SYMMLQ.
% SYMMLQ is an implementation of symmetric cg that applies to
% any symmetric A and will converge more rapidly than LSQR.
% If A is positive definite, there are other implementations of
% symmetric cg that require slightly less work per iteration
% than SYMMLQ (but will take the same number of iterations).
%
% Notation
% --------
%
% The following quantities are used in discussing the parameters:
%
% Abar   =  (   A    ),          bbar  =  ( b )
%           ( damp*I )                    ( 0 )
%
% r      =  b  -  A*x,           rbar  =  bbar  -  Abar*x
%
% rNorm  =  sqrt( norm(r)**2  +  damp**2 * norm(x)**2 )
%        =  norm( rbar )
%
% relpr  =  the relative precision of floating-point arithmetic
%           on the machine being used.  On most machines,
%           relpr is about 1.0d-16 in double precision.
%
% Tomlab LSQR (Tlsqr)  minimizes the function rNorm with respect to x.
%
% function [ x, iStop, Iter, rNorm, xNorm, StdErr, aNorm, aCond, arNorm ] =  ...
%     Tlsqr( m, n, Aname, iw, rw, b, damp, aTol, bTol, condLim, MaxIter, ...
%            WantStdErr, PrintFile, D, x_0 )
%
% INPUT: Tlsqr needs at least six parameters. The default values given below
%        are set by the MEX-interface before calling the LSQR solver
%
%  m          Number of rows in A
%
%     ------  ADVANCED MEMORY HANDLING: --------
%    m(2)     Optionally another element may be added to the parameter m input
%             to initiate advanced memory handling.
%             This is applicable when doing many calls to Tlsqr with the
%             same memory demand each time.
%
%             When m has length 2, 2nd parameter is Alloc (default 0)
%
%             Alloc 0  - Allocate and deallocate memory as usual
%             Alloc 1  - Only allocate memory - no return parameters
%             Alloc 2  - Allocate memory, solve problem, no deallocate
%             Alloc 3  - Solve problem, no allocate or deallocate
%             Alloc 4  - Solve problem, deallocate
%             Alloc 5  - Only deallocate memory - no return parameters
%
%             Standard loop of N calls:
%             Normal sequence of calls if x_0 is never (or always) set
%             Do 1st call in loop with Alloc=2
%             Do 2nd call in loop with Alloc=3
%              ...
%             Do N-1th call in loop with Alloc=3
%             Do Nth call in loop with Alloc=4
%
%             If last call in loop was with Alloc=3, instead:
%             Do final call with Alloc=5, no output
%
%             If Warm start is to be used for the 2nd, 3rd etc. call,
%             remember to either set x_0 = zeros(n,1) as input for the
%             1st call or do the following:
%
%             Alternative warm start (loop of N calls):
%             Do one initial call with Alloc=1, x_0 zero vector, no output
%             Do 1st call in loop with Alloc=3, x_0 empty
%             Do 2nd call in loop with Alloc=3, x_0 set
%              ...
%             Do N-1th call in loop with Alloc=3, x_0 set
%             Do Nth call in loop with Alloc=4, x_0 set
%
%  n          Number of columns in A
%  Aname      The m by n matrix A, or a name of a callback routine
%             that computes A*x (mode = 1) and A'*x (mode = 2)
%             If Aname is a numerical matrix, all computations are done in MEX
%             As callback routine the Tomlab standard callback routine
%             Tlsqrmat.m may be used. Tlsqrmat either computes A*x and
%             A'*x in Matlab using A as iw (see 2. for iw)
%             or calls the user routine given as the string iw (see 3. for iw)
%  iw         1. Empty (not used) if Aname is a numerical matrix
%             2. Matrix A if Aname = 'Tlsqrmat' and the computation of the
%             matrix products should be done in Matlab
%             3. Name of a user function to be called from Tlsqrmat
%             4. Arbitrary user information sent as callback, if Aname
%             is a user function (not Tlsqrmat).
%
%  rw         1. Empty (not used) if Aname is a numerical matrix
%             2. Additional user information, sent to the user routine as last
%             argument if nonempty and Aname = 'tlsqrmat'.
%             3. Arbitrary user information sent as callback, if Aname
%             is a user function (not Tlsqrmat).
%
%  b          The m by 1 vector b
%
%  damp       The damping parameter for problem 3 above. (Default 0)
%             If the system A*x = b is incompatible, values of damp in the
%             range 0 to sqrt(relpr)*norm(A) will probably have a negligible
%             effect. Larger values of damp will tend to decrease the norm of x
%             and reduce the number of iterations required by LSQR.
%
%             The work per iteration and the storage needed by LSQR are the
%             same for all values of damp.
%
% aTol        An estimate of the relative error in the data defining the matrix
%             A. For example, if A is accurate to about 6 digits, set
%             aTol = 1.0e-6 .  Default value: 0
%
% bTol        An estimate of the relative error in the data defining the rhs
%             vector b.  For example, if b is accurate to about 6 digits, set
%             bTol = 1.0e-6 . Default value: 0
%
% condLim     An upper limit on cond(Abar), the apparent condition number of
%             the matrix Abar. Iterations will be terminated if a computed
%             estimate of cond(Abar) exceeds condLim. Intended to prevent
%             certain small or zero singular values of A or Abar from coming
%             into effect and causing unwanted growth in the computed solution.
%
%             condLim and damp may be used separately or
%             together to regularize ill-conditioned systems.
%
%             Normally, condLim should be in the range 1000 to 1/relpr.
%             Suggested value (Default 0):
%               condLim = 1/(100*relpr)  for compatible systems,
%               condLim = 1/(10*sqrt(relpr)) for least squares.
%
%             Note:  If the user is not concerned about the parameters
%             aTol, bTol and condLim, any or all of them may be set
%             to zero.  The effect will be the same as the values
%             relpr, relpr and 1/relpr respectively.
%
% MaxIter     An upper limit on the number of iterations.  Suggested value:
%             MaxIter = n/2 for well-conditioned systems
%                           with clustered singular values,
%             MaxIter = 4*n   otherwise.  Default = max (50,4*n)
%
% WantStdErr  An 0/1 flag to say if the array StdErr(*) of standard error
%             estimates should be computed.  If m > n  or  damp >  0,  the
%             system is overdetermined and the standard errors may be useful.
%             (See the first LSQR reference.).  Default 0 or if m == n, 0.
%
% PrintFile   Filename for Tlsqr printfile. Default if empty is no printfile.
%
% D           Preconditioning. Only n diagonal elements of D are given
%
% x_0         Initial solution vector. Normally LSQR starts with x=0, but
%             if x_0 is given (n-vector), LSQR tries a warm start with x_0
%
% OUTPUT:
% x(1:n)      Returns the computed solution x.
%
% StdErr      If WantStdErr is true, the dimension of StdErr must be n or more.
%             StdErr then returns standard error estimates for the components
%             of x. For each i, StdErr(i) is set to the value
%                   rNorm * sqrt( sigma(i,i) / t ),
%             where sigma(i,i) is an estimate of the i-th diagonal of the
%             inverse of Abar(transpose)*Abar
%                        and  t = 1      if  m .le. n,
%                             t = m - n  if  m .gt. n  and  damp = 0,
%                             t = m      if  damp .ne. 0.
%
%             If WantStdErr is false, StdErr is empty
%
% iStop       An integer giving the reason for termination:
%
%             0  x = 0  is the exact solution.  No iterations were performed.
%
%             1  The equations A*x = b are probably compatible.  Norm(A*x - b)
%                is sufficiently small, given the values of aTol and bTol.
%
%             2  damp is zero. The system A*x = b is probably not compatible.
%                A least-squares solution has been obtained that is
%                sufficiently accurate,  given the value of aTol.
%
%             3  damp is nonzero. A damped least-squares solution has been
%                obtained that is sufficiently accurate, given the value of aTol
%
%             4  An estimate of cond(Abar) has exceeded condLim. The system
%                A*x = b appears to be ill-conditioned.  Otherwise, there
%                could be an error in the internal matrix product routine.
%
%             5  The iteration limit MaxIter was reached.
%
% Iter        The number of iterations performed.
%
% aNorm       An estimate of the Frobenius norm of Abar. This is the
%             square-root of the sum of squares of the elements of Abar. If
%             damp is small and if the columns of A have all been scaled to
%             have length 1.0, aNorm should increase to roughly sqrt(n). A
%             radically different value for aNorm may indicate an error in
%             the matrix product routine (inconsistency between modes 1 and 2?).
%
% aCond       An estimate of cond(Abar), the condition number of Abar. A very
%             high value of aCond may indicate an error in the matrix product.
%
% rNorm       An estimate of the final value of norm(rbar), the function being
%             minimized (see notation above).  Small if A*x = b has a solution.
%
% arNorm      An estimate of the final value of norm( Abar(transpose)*rbar ),
%             the norm of the residual for the usual normal equations. This
%             should be small in all cases. (arNorm will often be smaller than
%             the true value computed from the output vector x.)
%
% xNorm       An estimate of the norm of the final solution vector x.
%
%
% LSQR uses an iterative (conjugate-gradient-like) method.
%
% References
% ----------
%
% C.C. Paige and M.A. Saunders,  LSQR: An algorithm for sparse
%      linear equations and sparse least squares,
%      ACM Transactions on Mathematical Software 8, 1 (March 1982), pp. 43-71.
%
% C.C. Paige and M.A. Saunders,  Algorithm 583, LSQR: Sparse
%      linear equations and least-squares problems,
%      ACM Transactions on Mathematical Software 8, 2 (June 1982), pp. 195-209.
%
% C.L. Lawson, R.J. Hanson, D.R. Kincaid and F.T. Krogh,
%      Basic linear algebra subprograms for Fortran usage,
%      ACM Transactions on Mathematical Software 5, 3 (Sept 1979), pp. 308-325
% ----------------------------------------------------------------------------
%
%
%     LSQR development:
%     22 Feb 1982: LSQR sent to ACM TOMS to become Algorithm 583.
%     15 Sep 1985: Final F66 version.  LSQR sent to "misc" in netlib.
%     13 Oct 1987: Bug (Robert Davies, DSIR).  Have to delete
%                     if ( (one + dabs(t)) .le. one ) GO TO 200
%                  from loop 200.  The test was an attempt to reduce
%                  underflows, but caused w(i) not to be updated.
%     17 Mar 1989: First F77 version.
%     04 May 1989: Bug (David Gay, AT&T).  When the second beta is zero,
%                  rNorm = 0 and
%                  test2 = arNorm / (aNorm * rNorm) overflows.
%                  Fixed by testing for rNorm = 0.
%     05 May 1989: Sent to "misc" in netlib.
%     14 Mar 1990: Bug (John Tomlin via IBM OSL testing).
%                  Setting rhbar2 = rhobar**2 + dampsq can give zero
%                  if rhobar underflows and damp = 0.
%                  Fixed by testing for damp = 0 specially.
%     15 Mar 1990: Converted to lower case.
%     21 Mar 1990: d2norm introduced to avoid overflow in numerous
%                  items like  c = sqrt( a**2 + b**2 ).
%     04 Sep 1991: WantStdErr added as an argument to LSQR, to make
%                  standard errors optional.  This saves storage and
%                  time when StdErr(*) is not wanted.
%     13 Feb 1992: iStop now returns a value in [1,5], not [1,7].
%                  1, 2 or 3 means that x solves one of the problems
%                  Ax = b,  min norm(Ax - b)  or  damped least squares.
%                  4 means the limit on cond(A) was reached.
%                  5 means the limit on iterations was reached.
%     07 Dec 1994: Keep track of dxmax = max_k norm( phi_k * d_k ).
%                  So far, this is just printed at the end.
%                  A large value (relative to norm(x)) indicates
%                  significant cancellation in forming
%                  x  =  D*f  =  sum( phi_k * d_k ).
%                  A large column of D need NOT be serious if the
%                  corresponding phi_k is small.
%     27 Dec 1994: Include estimate of alfa_opt in iteration log.
%                  alfa_opt is the optimal scale factor for the
%                  residual in the "augmented system", as described by
%                  A. Bjorck (1992),
%                  Pivoting and stability in the augmented system method,
%                  in D. F. Griffiths and G. A. Watson (eds.),
%                  "Numerical Analysis 1991",
%                  Proceedings of the 14th Dundee Conference,
%                  Pitman Research Notes in Mathematics 260,
%                  Longman Scientific and Technical, Harlow, Essex, 1992.
%
%     Kenneth Holmstrom:
%     7  Oct 2000: LSQR called from Matlab using MEX-file interface
%     19 Jan 2003: Revised for Tomlab v4.0, including different callbacks
%                  Added warm start
%
%     Michael A. Saunders                  mike@sol-michael.stanford.edu
%     Dept of Operations Research          na.Msaunders@na-net.ornl.gov
%     Stanford University
%     Stanford, CA 94305-4022              (415) 723-1875
%
% function [ x, iStop, Iter, rNorm, xNorm, StdErr, aNorm, aCond, arNorm ] =  ...
%     Tlsqr( m, n, Aname, iw, rw, b, damp, aTol, bTol, condLim, MaxIter, ...
%            WantStdErr, PrintFile, x_0 )
%
%
% The above comments are a modified version of the original comments
% by Michael A. Saunders, suitably changed by Kenneth Holmstrom
% for the Matlab version of LSQR
%-----------------------------------------------------------------------

%# mex

function [ x, iStop, Iter, rNorm, xNorm, StdErr, aNorm, aCond, arNorm ] =  ...
    Tlsqr( m, n, Aname, iw, rw, b, damp, aTol, bTol, condLim, MaxIter, ...
           WantStdErr, PrintFile, D, x_0 )

help Tlsqr;