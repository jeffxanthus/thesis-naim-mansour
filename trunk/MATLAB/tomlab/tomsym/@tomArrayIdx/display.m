function display(p)
% tomArrayIdx/display - Command window display of a tomArrayIdx

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

disp(' ');
disp([inputname(1),' = tomSymIdx:'])
disp(' ');
disp(['   ' char(p) ' in ' mat2str(double(p)) ]);
disp(' ');
