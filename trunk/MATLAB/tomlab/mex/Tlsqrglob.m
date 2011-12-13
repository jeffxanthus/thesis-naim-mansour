% function y = Tlsqrglob( mode, m, n, x, Aname, rw )
%
% Callback routine for Tlsqr, for use with global matrices.
% The user must copy and edit this function to suit the 
% application. Sample code is given below.
%
% INPUT:
% mode   =1 Return y = A*x, =2 Return y = A' * x
% m         Number of rows in A
% n         Number of columns in A
% x         Vector of length n (mode=1)
%                         or m (mode=2), to be multiplied
% Aname     Not used
% rw        Not used
%
%
% Example of use:
%
% global A     % Declare before setting A
%
% A = [ coefficient matrix ]
% b = [ rhs ]
%
% iw=[]; rw=[];
% 
% [x, ...] = Tlsqr(m, n, 'Tlsqrglob', iw, rw, b, ...)
% 
% This will make Tlsqr call Tlsqrglob for each multiplication:
%
%   A*x   (if mode==1)
%   A'*x  (if mode==2)

function y = Tlsqrglob( mode, m, n, x, Aname, rw )

% IMPORTANT: 
% Declare global matrix/matrices here as needed for
% the calculation of A*x

global A

% The user must, if necessary, modify the code below to 
% perform matrix*vector multiplication along the lines of:

if mode==1
 y = A*x;
else
 y = A'*x;
end