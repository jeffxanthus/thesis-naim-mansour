% function f = funch2str(f)
%
% Converts a function handle to a string
%
% INPUTS:
%         
% f     Variable that may be a function handle
%
% OUTPUTS:
%
% f     The string representing the function handle

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.1.0$
% Written Dec 16, 2005.   Last modified Dec 16, 2005.

function f = funch2str(f)

if isa(f,'function_handle')
   f_func = functions(f);
   f = f_func.function;
end

% MODIFICATION LOG
%
% 051216  med  Written