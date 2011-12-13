function c = mcptest1_c(x, Prob)

% initialize
x = x(:);

% obtain variable values
z = x(1:6);     % quantities
p_d = x(7:9);   % demand price
p_s = x(10:11); % supply price

% check for domain violation
if (p_d <= 0)
  domerr = 1;
  return;
end

cost   = [ 0.225; 0.153; 0.162; 0.225; 0.162; 0.126; ];
demand = [ -325; -300; -275; ] ./ (p_d.^0.5);
supply = [ 325; 575; ];

A = sparse( [
   1   0   0   1   0   0;
   0   1   0   0   1   0;
   0   0   1   0   0   1;
  -1  -1  -1   0   0   0;
   0   0   0  -1  -1  -1;
] );

c = [ (-A'*[p_d; p_s] + cost); A*z + [demand; supply] ];