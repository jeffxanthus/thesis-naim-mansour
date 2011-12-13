% Operator precedence used to determine when to introduce
% parenthesis in displayed code. This should be the same as the Matlab
% operator precedence.
%
% 0 - Constants and variables
% 1 - Transpose (.'), power (.^), complex conjugate transpose ('), matrix power (^)
% 2 - Unary plus (+), unary minus (-), logical negation (~)
% 3 - Multiplication (.*), right division (./), left division (.\), matrix
%     multiplication (*), matrix right division (/), matrix left division (\)
% 4 - Addition (+), subtraction (-)
% 5 - Colon operator (:)
% 6 - Less than (<), less than or equal to (<=), greater than (>), greater than
%     or equal to (>=), equal to (==), not equal to (~=)
% 7 - Element-wise AND (&)
% 8 - Element-wise OR (|)
