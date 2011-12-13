% @TOMSYM
%
% Files
%   abs                       - tomSym/abs - Overloaded function 
%   acos                      - tomSym/acos - Overloaded function
%   acosh                     - tomSym/acosh - Overloaded function
%   acot                      - tomSym/acot - Overloaded function
%   acoth                     - tomSym/acoth - Overloaded function
%   acsc                      - tomSym/acsc - Overloaded function
%   acsch                     - tomSym/acsch - Overloaded function
%   angle                     - tomSym/angle - Overloaded function
%   asec                      - tomSym/asec - Overloaded function
%   asech                     - tomSym/asech - Overloaded function
%   asin                      - tomSym/asin - Overloaded function
%   asinh                     - tomSym/asinh - Overloaded function
%   atan                      - tomSym/atan - Overloaded function
%   atan2                     - tomSym/atan2 - Overloaded function
%   atanh                     - tomSym/atanh - Overloaded function
%   bsxfun                    - tomSym/bsxfun - Overloaded function
%   cat                       - tomSym/cat - Overloaded function
%   ceil                      - tomSym/ceil - Overloaded function
%   char                      - tomSym/char - Convert tomSym to char.
%   complementary             - complementary - overloaded function
%   conj                      - tomSym/conj - Overloaded function 
%   constpart                 - tomSym/constpart - The zeroth term of the McLaurin series of a tomSym
%   cos                       - tomSym/cos - Overloaded function
%   cosh                      - tomSym/cosh - Overloaded function
%   cot                       - tomSym/cot - Overloaded function
%   csc                       - tomSym/csc - Overloaded function
%   ctranspose                - tomSym/ctranspose - Overloaded operator
%   cumsum                    - tomSym/cumsum - Overloaded function
%   dblquad                   - tomSym/dblquad - Numeric double quadrature
%   derivative                - tomSym/derivative - The symbolic derivative of a tomSym object.
%   derivatives               - tomSym/derivatives - The symbolic derivatives of a tomSym object
%   det                       - tomSym/det - Overloaded function
%   diag                      - tomSym/diag - Overloaded function
%   diff                      - tomSym/diff - Difference - Overloaded function
%   DiracDelta                - tomSym/DiracDelta - Placeholder derivative of a discontinuous function.
%   DiracDeltas               - tomSym/DiracDeltas - Placeholder derivative of a function with multiple discontinuities.
%   disp                      - tomSym/disp - Command window display of a tomSym
%   display                   - tomSym/display - Command window display of a tomSym
%   double                    - tomSym/double - Convert tomSym constant to double
%   end                       - tomSym/end - Overloaded - Get the last index of a tomSym object
%   eps                       - tomSym/rem - Overload the remainder function
%   eq                        - tomSym/eq - Overloaded == operator
%   erf                       - tomSym/erf - Overloaded function
%   erfc                      - tomSym/erfc - Overloaded function
%   erfcinv                   - tomSym/erfcinv - Overloaded function
%   erfinv                    - tomSym/erfinv - Overloaded function
%   estpattern                - tomSym/estpattern - Estimate the sparsity pattern of a tomSym
%   eval                      - tomSym/eval - Evaluate a tomSym object in the current workspace
%   exp                       - tomSym/exp - Overloaded function
%   extractConstraints        - extrectContraitns - Move constraints from subjectTo into constraint set.
%   eye                       - tomSym/eye - Overloaded function
%   ezplot                    - tomSym/ezplot - Function plot
%   feval                     - tomSym/feval - tomSym-compatible function call.
%   fix                       - tomSym/fix - Overloaded function
%   fliplr                    - tomSym/fliplr - Overloaded function
%   flipud                    - tomSym/flipud - Overloaded function
%   floor                     - tomSym/floor - Overloaded function
%   full                      - tomSym/full - Overloaded function
%   fzero                     - tomSym/fzero - same as tomfzero
%   gamma                     - tomSym/gamma - Overloaded function
%   gammainc                  - tomSym/gammainc - Overloaded function
%   gammaln                   - tomSym/gammaln - Overloaded function
%   ge                        - tomSym/ge - Overloaded >= operator
%   getdiag                   - tomSym/getdiag - Getdiag for tomSym.
%   gt                        - tomSym/gt - Overloaded > operator
%   holdInterpolationMatrix   - holdInterpolationMatrix - Value matrix for interpPoly
%   horzcat                   - tomSym/horzcat - Overload the [] operator
%   hownonlinear              - hownonlinear - overloaded function
%   ifThenElse                - tomSym/ifThenElse - Symbolic if/then/else
%   imag                      - tomSym/imag - Overloaded function
%   interp1                   - tomSym/interp1 - Overloaded function
%   interp1h                  - tomSym/interp1h - Overloaded function
%   interp1l                  - tomSym/interp1l - Overloaded function
%   interp1s                  - tomSym/interp1s - Overloaded function
%   interp1sDot               - tomSym/interp1sDot - Derivative of interp1s, overloaded
%   interp2                   - TOMSYM/INTERP2 - overloaded function
%   interpn                   - TOMSYM/INTERPN - overloaded function
%   inv                       - tomSym/inv - Overloaded function
%   isdependent               - tomSym/isdependent - Determine if f is dependent of the symbol x
%   isempty                   - tomSym/isempty - Overloaded function
%   isequal                   - tomSym/isequal - Overloaded function
%   isfinite                  - tomSym/isfinite - Overloaded function
%   isinf                     - tomSym/isinf - Overloaded function
%   isnan                     - tomSym/isnan - Overloaded function
%   isreal                    - tomSym/isreal - Overloaded function
%   istomsymbol               - tomSym/istomsymbol - Returns "true" if a is a tomSym symbol.
%   kkt1                      - kkt1 - First order Karush-Kuhn-Tucker conditions
%   kron                      - tomSym/kron - Overload the Kronecker product
%   ldivide                   - tomSym/ldivide - Overload the .\ left division operator
%   le                        - tomSym/ge - Overloaded <= operator
%   length                    - tomSym/length - Get the length of a tomSym object
%   linearInterpolationMatrix - linearInterpolationMatrix - Value matrix for interpPoly
%   linspace                  - tomSym/linspace - Overloaded function
%   log                       - tomSym/log - Overloaded function
%   log10                     - tomSym/log10 - Overloaded function
%   log2                      - tomSym/log2 - Overloaded function
%   logm                      - tomSym/logm - Overloaded function
%   lookup                    - tomSym/lookup - Overloaded function
%   lt                        - tomSym/lt - Overloaded < operator
%   mat2str                   - tomSym/mat2str - Convert tomSym to a string.
%   max                       - tomSym/max - Overloaded
%   maxIndicator              - maxIndicator - Overloaded function
%   mcode                     - tomSym/mcode - Generate m-code from a tomSym object.
%   mcodestr                  - tomSym/mcodestr - Convert tomSym to m-code
%   MI                        - tomSym/MI - Matrix Inequality
%   min                       - tomSym/min - Overloaded
%   minus                     - tomSym/minus - Overload the minus operator
%   mldivide                  - tomSym/mrdivide - Overload the left divide operator
%   mod                       - tomSym/mod - Overloaded function
%   monofun                   - monofun - Overloaded function
%   monoinv                   - monoinv - Overloaded function
%   mpower                    - tomSym/mpower - Overload the mpower operator
%   mrdivide                  - tomSym/mrdivide - Overload the divide operator
%   mtimes                    - tomSym/mtimes - Overload the multiplication operator
%   nnz                       - tomSym/spy - Number of nonzero matrix elements
%   nOperands                 - tomSym/nOperands - Get number of operands from a tomSym
%   norm                      - tomSym/norm - Overloaded function
%   not                       - tomSym/not - Overloaded function
%   num2str                   - tomSym/num2str - Convert tomSym to a string.
%   numel                     - tomSym/numel - Returns prod(size(a)) for tomSym a.
%   operand                   - tomSym/operand - Get an operand from a tomSym
%   operands                  - tomSym/operands - Get all operands from a tomSym
%   operator                  - tomSym/operator - Get the operator from a tomSym
%   pattern                   - tomSym/pattern - The sparsity pattern of a tomSym object
%   plus                      - tomSym/plus - Overload the plus operator
%   positiveSemidefinite      - tomSym/positiveSemidefinite - Overloaded function
%   power                     - tomSym/power - Overload the mpower operator
%   ppnval                    - TOMSYM/PPNVAL - Overloaded function
%   ppval                     - tomSym/ppval - Overloaded function
%   prod                      - tomSym/prod - Overloaded function 
%   prodJ1                    - tomSym/prodJ1 - overloaded function 
%   prodJ1J1                  - tomSym/prodJ1J1 - overloaded function 
%   psi                       - tomSym/psi - Overloaded function
%   quad                      - tomSym/quad - Numeric quadrature
%   rdivide                   - tomSym/rdivide - Overload the ./ division operator
%   real                      - tomSym/real - Overloaded function
%   rem                       - tomSym/rem - Overload the remainder function
%   repmat                    - tomSym/repmat - Overloaded function
%   reshape                   - tomSym/reshape - Overloaded function
%   rewriteV                  - rewriteV - Rewrite an optimization problem to avoid sharp corners.
%   round                     - tomSym/round - Overloaded function
%   scalecolumns              - tomSym/scalecolumns - Overloaded function
%   scalerows                 - tomSym/scalerows - Overloaded function
%   sec                       - tomSym/sec - Overloaded function
%   setdiag                   - tomSym/setdiag - Setdiag for tomSym.
%   setSymmetric              - tomSym/setSymmetric - Create a symmetric matrix
%   sign                      - tomSym/sign - Overloaded function 
%   sin                       - tomSym/sin - Overloaded function
%   sinh                      - tomSym/sinh - Overloaded function
%   size                      - tomSym/size - Get the size of a tomSym object
%   smplus                    - tomSym/smplus - Scalar + matrix addition
%   smtimes                   - tomSym/smtimes - Scalar x matrix multiplication
%   sparse                    - tomSym/sparse - Overloaded function
%   spdiags                   - tomSym/spdiags - Overloaded function
%   spower                    - tomSym/spower - Element-wise power to a scalar exponent.
%   spy                       - tomSym/spy - Visualize sparsity pattern for a tomSym
%   sqrt                      - tomSym/sqrt - Overloaded function
%   sub2ind                   - tomSym/sub2ind - Overloaded function.
%   submatrix                 - tomSym/submatrix - Overloaded function
%   subsasgn                  - tomSym/subsasgn - Subscripted assignemnt, overloaded for tomSym.
%   subsindex                 - tomSym/subsindex - Not possible
%   subsref                   - tomSym/subsref - Object index lookup
%   subststruct               - tomSym/subststruct - Substitue tomSym symbols from a struct
%   subsymb                   - tomSym/subsymb - Get a sub-symbol from a tomSym
%   sum                       - tomSym/sum - Overloaded function
%   sym                       - tomSym/sym - convert a tomSym object to a symbolic toolbox object
%   sym2eig                   - sym2eig Convert symbolic constraints into an eigenvalue problem
%   sym2prob                  - sym2prob - Compile symbolic function/constraints into a Prob struct.
%   symbols                   - tomSym/symbols - List all symbols used in a tomSym
%   tan                       - tomSym/tan - Overloaded function
%   tanh                      - tomSym/tanh - Overloaded function
%   times                     - tomSym/times - Overload the .* multiplication operator (Hadamard product)
%   tomCmp                    - tomCmp - Test if the operator of a tomSym is of a certain type
%   tomnumeric                - tomSym/tomnumeric - isnumeric for tomSym
%   tomSym                    - tomSym/tomSym - Class constructor
%   tomUnGlobalize            - tomSym/tomUnGlobalize - undo the effect of tomGlobalize
%   trace                     - tomSym/trace - Overloaded function
%   transpose                 - tomSym/transpose - Overloaded operator
%   tril                      - tomSym/tril - Overloaded function 
%   triu                      - tomSym/triu - Overloaded function 
%   uminus                    - tomSym/uminus - Overloaded function
%   uplus                     - tomSym/uplus - Overloaded function
%   vec                       - tomSym/vec - Transform a tomSym object into a column vector
%   vertcat                   - tomSym/vertcat - Overload the [;] operator
%   wrap                      - tomSym/wrap - Overloaded
%   wrapJ                     - tomSym/wrapJ - Overloaded
