function testfminunc

fprintf('\n');
fprintf('------------------NEW TEST ----------------------------------\n');
fprintf('\n');
disp('Test of fminunc, from OPT TB 2.0 page 4.62. Simple f=sin(x)+3');
fprintf('\n');

x_0 = [2];

% Gradient is supplied
if exist('optimset')
   options=optimset;
else
   options=[];
end
options.GradObj='on';
options.LargeScale='off';

fprintf('Test problem in file %s','fminunc_fg');
fprintf('\n');   
fprintf('\n');   
fprintf('Gradient is supplied. Starting point x = 2.\n');
fprintf('\n');   

options.Display='iter';

[x,f_k,ExitFlag,Output] = fminunc('fminunc_fg',x_0,options);

format compact
disp('Output Structure')
disp(Output)
xprinte(x,'x:             ');
fprintf('\n');   
fprintf('Function value %30.20f. ExitFlag %d. ',f_k,ExitFlag);
fprintf('Iterations %d.\n',Output.iterations);

fprintf('\n');
fprintf('------------------NEW TEST ----------------------------------\n');
fprintf('\n');
disp('Test of fminunc. Example from OPT TB 2.0 page 4.61-4.62');

x_0 = [1;1];

% Gradient is supplied
if exist('optimset')
   options=optimset;
else
   options=[];
end
options.GradObj='on';
options.LargeScale='off';
fprintf('\n');   
fprintf('Test problem in file %s','fminunc_fg2');
fprintf('\n');   
fprintf('\n');   
fprintf('Gradient is supplied. ');
fprintf('Starting point x = (1,1).\n');
fprintf('\n');   

options.Display='iter';

[x,f_k,ExitFlag,Output] = fminunc('fminunc_fg2',x_0,options);

format compact
disp('Output Structure')
disp(Output)
xprinte(x,'x:             ');
fprintf('\n');   
fprintf('Function value %30.20f. ExitFlag %d. ',f_k,ExitFlag);
fprintf('Iterations %d.\n',Output.iterations);

fprintf('\n');
fprintf('------------------NEW TEST ----------------------------------\n');
fprintf('\n');

disp('Test of fminunc, from OPT TB 2.0 page 4.61-4.62');
fprintf('\n');   
fprintf('Test problem in file %s','fminunc_fg2');
fprintf('\n');   
fprintf('\n');   
fprintf('Use numerical differences for the gradient. ');
fprintf('Starting point x = (1,1).\n');
fprintf('\n');   

x_0 = [1;1];

% Gradient is NOT supplied (even if defined in file)
if exist('optimset')
   options=optimset;
else
   options=[];
end
options.GradObj='off';
options.LargeScale='off';
options.Display='iter';
options.Diagnostics='on';

[x,f_k,ExitFlag,Output] = fminunc('fminunc_fg2',x_0, options);

format compact
disp('Output Structure')
disp(Output)
xprinte(x,'x:             ');
fprintf('\n');   
fprintf('Function value %30.20f. ExitFlag %d. ',f_k,ExitFlag);
fprintf('Iterations %d.\n',Output.iterations);