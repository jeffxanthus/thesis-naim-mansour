% function y = Tlsqrmat( mode, m, n, x, Aname, rw )
%
% Callback routine for Tlsqr
%
% INPUT:
% mode   =1 Return y = A*x, =2 Return y = A' * x
% m         Number of rows in A
% n         Number of columns in A
% x         Vector of length n (mode=1) or m (mode=2), to be multiplied
% Aname     A matrix, or name of function that computes A*x, or A'*x
% rw        Extra parameter, not used
%
% If A is an explicit matrix and the product is to be computed in Matlab, do
%
% rw = []; [x, ...] = Tlsqr(m, n, 'Tlsqrmat', A, rw, b, ...)
% 
% then Tlsqrmat is doing A*x, or A'*x (using Aname as A).
%
% --------------------------------------------------------------------------
% Doing the call:
%
% rw = []; [x, ...] = Tlsqr(m, n, 'Tlsqrmat', 'UserFile', rw, b)
%
% Tlsqrmat will make the following call to UserFile:
%
% y = feval('UserFile', mode, m, n, x) 
%
% --------------------------------------------------------------------------
% Setting rw nonempty and doing the call:
%
% [x, ...] = Tlsqr(m, n, 'Tlsqrmat', 'UserFile', rw, b)
%
% will make Tlsqrmat add rw as extra parameter in the call to UserFile:
%
% y = feval('UserFile', mode, m, n, x, rw)

function y = Tlsqrmat( mode, m, n, x, Aname, rw )

if ischar(Aname)
   if isempty(rw)
      y = feval( Aname, mode, m, n, x );
   else
      y = feval( Aname, mode, m, n, x, rw );
   end
else
   if mode==1,  y = Aname*x;  else  y = Aname'*x;  end
end