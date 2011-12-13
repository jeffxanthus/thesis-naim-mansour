function testfgoalattain

fprintf('\n');
fprintf('------------------NEW TEST ----------------------------------\n');
fprintf('\n');
disp('Test of fgoalattain, from OPT TB 3');
fprintf('\n');

if exist('optimset')
   options=optimset;
else
   options=[];
end

options.Display = 'iter';

A =  [ -0.5  0  0;  0  -2  10;  0  1  -2 ];
B =  [ 1  0;  -2  2;  0  1 ];
C =  [ 1  0  0;  0  0  1 ];
goal = [-5, -3, -1];
weight = abs(goal);
x0 = [ -1 -1 -1 -1]; 
lb = repmat(-4,size(x0)) 
ub = repmat(4,size(x0))

eigfun = 'testfgoalattain_f';

[x,f_k,attainfactor,ExitFlag,output,lambda] = ...
        fgoalattain(eigfun,x0,goal,weight,[],[],[],[],lb,ub,[],options); 

format compact
fprintf('\n');   
xprinte(x,'x:             ');
fprintf('\n');   
fprintf('Function value %30.20f. ExitFlag %d.\n',f_k,ExitFlag);

options.GoalsExactAchieve = 3;

fprintf('\n');
fprintf('------------------NEW TEST ----------------------------------\n');
fprintf('\n');
disp('Test of fgoalattain, from OPT TB 3, Exact goals');
fprintf('\n');

% This problem has two minima, NPSOL finds the lowest one, the global minimum

[x,f_k,attainfactor,ExitFlag,output,lambda] = ...
    fgoalattain(eigfun,x0,goal,weight,[],[],[],[],lb,ub,[],options);

format compact
fprintf('\n');   
xprinte(x,'x:             ');
fprintf('\n');   
fprintf('Function value %30.20f. ExitFlag %d.\n',f_k,ExitFlag);
