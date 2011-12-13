% nsdp_x2xY.m:
%
% function [x,Y] = nsdp_x2xY(xY, NSDP)
%
% TOMLAB gateway routine for nonlinear semidefinite programming.
% nsdp_x2xY separates the standard variables from the matrix variables by
% using the information stored in the struct Prob.PENOPT.NSDP.
%
% The matrix variables Y are returned as a cell array.
%
% Bjorn Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.1.0$
% Written October 8, 2008. Last modified October 8, 2008.
function [x,Y] = nsdp_x2xY(xY, NSDP)

% Count standard and matrix variables
n_X = length(xY);
n_matrix_vars = length(NSDP);
for i = 1:n_matrix_vars
   n_X = n_X - length(NSDP(i).Y_var_index);
end

% Separate standard variables
x = xY(1:n_X);

% Matrix variables are stored after the standard ones.
Y_vector = xY(n_X+1:end);

% Step through the struct array NSDP and create
% a cell array of Y-matrices with updated values.
pos = 1;
for i = 1:n_matrix_vars
   % Indices of variable elements
   Yi = NSDP(i).Y_var_index;
   
   % Number of variable elements
   nYi = length(Yi);
   
   % Extract only upper triangle
   Y_i = triu(NSDP(i).Y);
   
   % Update upper triangle
   Y_i(Yi) = Y_vector(pos:pos+nYi-1);
   
   % Make Y symmetric
   Y{i} = Y_i + Y_i' - diag(diag(Y_i));
   
   pos = pos + nYi;
end

% MODIFICATION LOG
%
% 081008 bjo  Wrote file