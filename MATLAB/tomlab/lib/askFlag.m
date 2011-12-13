% Code sequence for switch flag ask. Used by xxxOpt routines.
%
% Assumes the variable ask to be defined, and reverses this flag.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written May 30, 1998.  Last modified June 21, 1999.

fprintf('\nSwitch: toggles ASK/NOT ASK in problem definitions\n\n');

if ask > 0
   ask=askInit;
else
   ask=1;
end
if ask==1
   fprintf('\nSwitch is now on: Ask questions in problem definition\n');
else
   fprintf('\nSwitch is now off: ');
   fprintf('Do NOT ask questions in problem definitions\n\n');
end