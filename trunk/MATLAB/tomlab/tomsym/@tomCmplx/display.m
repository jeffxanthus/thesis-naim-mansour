function display(p)
% tomCmplx/display - Command window display of a tomCmplx

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-12-23 by rutquist for TOMLAB release 7.7

disp(' ');
disp([inputname(1),' = tomCmplx(' num2str(size(p,1)) 'x' num2str(size(p,2)) '):'])
disp('Real part:');
disp(p.re);
disp('Imaginary part:');
disp(p.im);
