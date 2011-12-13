function dc = mcptest1_dc(x, Prob)

% initialize
x = x(:);
p_d = x(7:9);   % demand price

% check for domain violation
if (p_d <= 0)
  domerr = 1;
  return;
end

A = sparse( [
   1   0   0   1   0   0;
   0   1   0   0   1   0;
   0   0   1   0   0   1;
  -1  -1  -1   0   0   0;
   0   0   0  -1  -1  -1;
] );

dc = [ sparse(6,6) -A'; 
        A diag( [ -0.5 * [-325; -300; -275;] ./ (p_d.^1.5); 0; 0] ) ];