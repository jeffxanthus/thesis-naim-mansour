% example_ann - tomSym ANN demonstration
%
% In this example it is shown how tomSym can be used to train a simple
% artificial neural network that solves the "xor" problem.
%
% This is only an example - neural networks are not usually trained this
% way, but the example shows how easily backpropagation is implemented in
% tomSym.

Ninputs  = 2;
Nhidden  = 4;
Noutputs = 1;

W1 = tom('W1',Nhidden,Ninputs+1);  % Layer 1 weights and offset
W2 = tom('W2',Noutputs,Nhidden+1); % Layer 2 weights and offset

% Loop over all indata
objective = 0;
for i = [-1 1]
    for j = [-1 1]
        % Simulate the ANN.
        indata     = [i; j];
        hiddendata = tanh(W1*[indata; -1]);
        outdata    = tanh(W2*[hiddendata; -1]);
        answer     = i*j;
        err        = (outdata-answer)^2;
        objective  = objective + err;
    end
end

constraints = { -5 <= W1 <= 5, -5 <= W2 <= 5};

% Compile and solve problem
options = struct;
options.name = 'ANN';
options.solver = 'multiMin';
options.xInit = 20;
solution = ezsolve(objective,constraints,[],options);