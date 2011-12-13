% optim_fgH.m
%
% function f = optim_fgH(x, Prob)
%
% optim_fgH is used to implement the OPT TB 2.0 interface
%
% The function f is returned
%
% If the gradient is computed, it is stored in the global variable
% NLP_g (and the corresponding x value in NLP_xg)
%
% If the Hessian is computed, it is stored in the global variable
% NLP_H (and the corresponding x value in NLP_xH)
%
% optim_fgH is called from the TOMLAB gateway function nlp_f.

% Kenneth Holmstrom, Tomlab Optimization AB, E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2006 by Tomlab Optimization AB, Sweden. $Release: 5.2.0$
% Written July 29, 1999.     Last modified Feb 19, 2006.

function f = optim_fgH(xNew, Prob)

global NLP_xg NLP_g NLP_xH NLP_H

Func    = Prob.OPTTB.f;
aF      = Prob.OPTTB.funArgIn;
oF      = Prob.OPTTB.funArgOut;
NumDiff = Prob.NumDiff;
x       = Prob.OPTTB.x;
x(:)    = xNew;                   % To handle cases when x is a matrix
n       = length(x);

if Prob.OPTTB.M7 == 0

   if NumDiff >0 | Prob.ADObj > 0 
      % Gradient should not be computed
      if oF < 0
         f=eval(Func);
      else
         if aF > 1
            f=feval(Func, x, Prob.varargin{:} );
         else
            f=feval(Func, x);
         end
      end
      NLP_xg=[];
      NLP_g=[];
      NLP_xH=[];
      NLP_H=[];
   elseif NumDiff | Prob.ADObj < 0       
      % Hessian should not be computed
      if oF < 0
         f=eval(Func);
         NLP_g = eval(Prob.OPTTB.g);
      elseif oF > 1
         if aF > 1
            [f,NLP_g] = feval(Func, x, Prob.varargin{:});
         else
            [f,NLP_g] = feval(Func, x);
         end
         % Skip NLP_g if not of correct size
         if length(NLP_g)~= n, NLP_g = []; end  
      else
         if aF > 1
            f = feval(Func, x, Prob.varargin{:});
         else
            f = feval(Func, x);
         end
         NLP_g=[];
      end
      if isempty(NLP_g)
         NLP_xg=[];
      else
         NLP_xg=x(:);
      end
      NLP_xH=[];
      NLP_H=[];
   else
      if oF < 0
         f=eval(Func);
         NLP_g = eval(Prob.OPTTB.g);
         if oF == -2
            NLP_H = eval(Prob.OPTTB.H);
         else
            NLP_H = [];
         end
      else
         if oF > 2
            if aF > 1
               [f,NLP_g,NLP_H] = feval(Func, x, Prob.varargin{:});
            else
               [f,NLP_g,NLP_H] = feval(Func, x);
            end
            % Skip NLP_g if not of correct size
            if length(NLP_g) ~= n, NLP_g = []; end  
            % Skip NLP_H if not of correct size
            if size(NLP_H,1) ~= n | size(NLP_H,2)~= n, NLP_H = []; end  
         elseif oF > 1
            if aF > 1
               [f,NLP_g] = feval(Func, x, Prob.varargin{:});
            else
               [f,NLP_g] = feval(Func, x);
            end
            % Skip NLP_g if not of correct size
            if length(NLP_g) ~= n, NLP_g = []; end  
            NLP_H = [];
         else
            if aF > 1
               f = feval(Func, x, Prob.varargin{:});
            else
               f = feval(Func, x);
            end
            NLP_g = [];
            NLP_H = [];
         end
         HFunc = Prob.HessMult;
         if ~isempty(HFunc)
            n = Prob.N;
            y = zeros(n,1);
            Hinfo = NLP_H;
            NLP_H = sparse(n,n);
            % Sick way to obtain the full Hessian - works, but slow 
            for i=1:n
                y(i) = 1;
                NLP_H(:,i) = feval(HFunc, Hinfo, y, Prob.varargin{:});
                y(i) = 0;
            end
         end
      end 
      if isempty(NLP_g)
         NLP_xg=[];
      else
         NLP_xg=x(:);
      end
      if isempty(NLP_H)
         NLP_xH=[];
      else
         NLP_xH=x(:);
      end
   end

else
   if NumDiff >0 | Prob.ADObj > 0 
      % Gradient should not be computed
      if aF > 1
         f=Func( x, Prob.varargin{:} );
      else
         f=Func( x);
      end
      NLP_xg=[];
      NLP_g=[];
      NLP_xH=[];
      NLP_H=[];
   elseif NumDiff | Prob.ADObj < 0       
      % Hessian should not be computed
      if oF > 1
         if aF > 1
            [f,NLP_g]=Func( x, Prob.varargin{:} );
         elseif oF > 1
            [f,NLP_g]=Func( x);
         end
         % Skip NLP_g if not of correct size
         if length(NLP_g) ~= n, NLP_g = []; end  
      else
         if aF > 1
            f=Func( x, Prob.varargin{:} );
         elseif oF > 1
            f=Func( x);
         end
         NLP_g=[];
      end
      if isempty(NLP_g)
         NLP_xg=[];
      else
         NLP_xg=x(:);
      end
      NLP_xH=[];
      NLP_H=[];
   else
      if oF > 2
         if aF > 1
            [f,NLP_g,NLP_H] = Func( x, Prob.varargin{:} );
         else
            [f,NLP_g,NLP_H] = Func( x);
         end
         % Skip NLP_g if not of correct size
         if length(NLP_g) ~= n, NLP_g = []; end  
         % Skip NLP_H if not of correct size
         if size(NLP_H,1) ~= n | size(NLP_H,2)~= n, NLP_H = []; end  
      elseif oF > 1
         if aF > 1
            [f,NLP_g]=Func( x, Prob.varargin{:} );
         else
            [f,NLP_g]=Func( x);
         end
         % Skip NLP_g if not of correct size
         if length(NLP_g) ~= n, NLP_g = []; end  
         NLP_H = [];
      else
         if aF > 1
            f=Func( x, Prob.varargin{:} );
         else
            f=Func( x);
         end
         NLP_g = [];
         NLP_H = [];
      end
      HFunc = Prob.HessMult;
      if ~isempty(HFunc)
         n = Prob.N;
         y = zeros(n,1);
         Hinfo = NLP_H;
         NLP_H = sparse(n,n);
         % Sick way to obtain the full Hessian - works, but slow 
         for i=1:n
             y(i) = 1;
             NLP_H(:,i) = feval(HFunc, Hinfo, y, Prob.varargin{:});
             y(i) = 0;
         end
      end
      if isempty(NLP_g)
         NLP_xg=[];
      else
         NLP_xg=x(:);
      end
      if isempty(NLP_H)
         NLP_xH=[];
      else
         NLP_xH=x(:);
      end
   end
end

% MODIFICATION LOG:
%
% 030114 hkh Major revision, avoid check of nargin(Func) here
% 030115 hkh Add handling of cell, inline and function_handle input
% 030126 hkh Handle the case when x is a matrix, and not a vector
% 030128 hkh Handle the HessMult option, generate Hessian
% 031101 hkh Change AutoDiff to new field ADObj, add test for Hessian for MAD
% 041119 hkh Test oF, allow both two or three outputs from user function
% 041222 hkh Handle new type of function handle in Matlab 7.x
% 060218 hkh Check for only one output argument, oF==1, avoid crash
% 060218 hkh Check for only one output argument, oF==1, avoid crash
% 060219 hkh Check size of NLP_g and NLP_H returned, if bad output from users
