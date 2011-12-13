% tomSym/feval - tomSym-compatible function call.
%
% FEVAL(FUN,x1,...,xn) calls FUN(x1,...,xn) if FUN is a function that is
% compatible with tomSym.
%
% If FUN is not compatible with tomSym, then FEVAL calls it via the tomSym
% WRAP function. FEVAL then attempts to automatically deduce all the
% information that is required by WRAP by calling FUN using random indata.
%
% The symbolic arguments will be substituted with random numbers when
% making test calls to FUN. If this results in an error, then WRAP must be
% used instead of feval. Additionaly, because FEVAL tries to estimate the
% sparsity pattern of FUN and its derivatives from these trial calls, those
% patterns may be erroneous, which will result in incorrect numeric
% derivatives. If this is the case then WRAP must be used instead of feval.
%
% For derivatives, TomSym will first look for m-files named on the form
% "funJ1", then attempt to use MAD, and if that fails, fall back to numeric
% differentiation.
%
% [y1,..,yn] = FEVAL(FUN,x1,...,xn) can be used if fun returns multiple
% output arguments.
%
% See also: wrap
