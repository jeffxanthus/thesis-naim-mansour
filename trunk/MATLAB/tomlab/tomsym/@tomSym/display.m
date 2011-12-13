function display(p)
% tomSym/display - Command window display of a tomSym

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

disp(' ');
disp([inputname(1),' = tomSym(' num2str(size(p,1)) 'x' num2str(size(p,2)) '):'])
disp(' ');
disp(['   ' char(p)]);
disp(' ');
