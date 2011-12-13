% tomSym/mcode - Generate m-code from a tomSym object.
%
% [code, tempD, header] = mcode(f, fname) generates m-code, representing an
% object f.
%
% All constants are moved to a tempD object, which needs to be supplied in
% the function call.
%
% The returned function header lists all symbols alphabetically, and the
% data input last.
