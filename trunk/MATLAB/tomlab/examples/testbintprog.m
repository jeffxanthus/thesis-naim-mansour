function testbintprog

fprintf('\n');
fprintf('------------------NEW TEST ----------------------------------\n');
fprintf('\n');
disp('Test of bintprog, from OPT TB 3');
fprintf('\n');

if exist('optimset')
   options=optimset;
else
   options=[];
end

options.Display='iter';
options.Display='final';

f = -ones(40,1);
A = sparse(4,40);
for i=1:4
    A(i,[(i-1)*10+1:i*10]) = -ones(1,10);
end
b = -[2 2 2 2]';
Aeq = ones(1,40);
beq = 10;

[x,f_k,ExitFlag] = bintprog(f,A,b,Aeq,beq);

format compact
fprintf('\n');   
xprinti(x,'x:             ');
fprintf('\n');   
fprintf('Function value %30.20f. ExitFlag %d.\n',f_k,ExitFlag);