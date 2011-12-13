% tomFiles sets the names of the m-Files to be used to compute the
% functions needed for optimization.
% 
% function Prob = tomFiles(Prob, f, g, H, c, dc, d2c, r, J, d2r, fc, gdc);
%
% INPUT:
%   Prob   The TOMLAB Problem structure
%             Fields used:
%
%   f      Name of function to compute function value f
%   g      Name of function to compute gradient vector g
%   H      Name of function to compute Hessian matrix H
%   c      Name of function to compute constraint vector c
%   dc     Name of function to compute constraint Jacobian dc
%   d2c    Name of function to compute 2nd part of 2nd der of Lagrangian
%   r      Name of function to compute residual vector r
%   J      Name of function to compute Jacobian matrix J
%   d2r    Name of function to compute 2nd part of 2nd derivative of 0.5*r'*r
%   fc     Name of function to compute function value f and constraint vector c
%   gdc    Name of function to compute gradient vector g and 
%          constraint Jacobian dc
%
% OUTPUT
%   Prob   The TOMLAB Problem structure
%             Fields used:
%             FUNCS.f
%             FUNCS.g
%             FUNCS.H
%             FUNCS.c
%             FUNCS.dc
%             FUNCS.d2c
%             FUNCS.r
%             FUNCS.J
%             FUNCS.d2r
%             FUNCS.fc
%             FUNCS.gdc
%             FUNCS.rJ  is set as []

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.5.0$
% Written Oct 10, 1998.    Last modified Aug 17, 2006.

function Prob = tomFiles(Prob, f, g, H, c, dc, d2c, r, J, d2r, fc, gdc)

if nargin < 14
   cdc=[];
   if nargin < 13
      fg=[];
      if nargin < 12
         gdc=[];
         if nargin < 11
            fc=[];
            if nargin < 10
               d2r=[];
               if nargin < 9
                  J=[];
                  if nargin < 8
                     r=[];
                     if nargin < 7
                        d2c=[];
                        if nargin < 6
                           dc=[];
                           if nargin < 5
                              c=[];
                              if nargin < 4
                                 H=[];
                                 if nargin < 3
                                    g=[];
end, end, end, end, end, end, end, end, end, end, end, end 

Prob.FUNCS = struct('f',f, 'g',g, 'H',H, 'c',c, 'dc',dc, 'd2c',d2c, ...
                   'r',r, 'J',J, 'd2r', d2r, 'fc', fc, 'gdc', gdc, ...
                   'fg', fg,'cdc', cdc,'rJ',[]);

% Here we could, safeguard, check if the routines are in the path etc.


% MODIFICATION LOG:
%
% 981013  mbk  Changed J=[] to H=[] if nargin < 4.
% 981023  hkh  Added new variable d2r for 2nd derivative of LS problems
% 990630  hkh  Write mideva.m with routines to be called
% 990919  hkh  Avoid writing to mideva.m if not Mideva is used
% 030524  hkh  Add fields fc and gdc for simulation problems
% 040102  hkh  Put mideva part in comments
% 040728  med  Name changed to tomFiles
% 050117  med  mlint revision
% 050623  ango Add fg and cdc fields
% 060814  med  FUNCS used for callbacks instead
% 060817  hkh  Must also set field rJ empty
