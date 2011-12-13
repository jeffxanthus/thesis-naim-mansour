% askparam prompts the user for 1 real value z in the interval [zLow,zUpp]
% using either inputR or inputR2 (which displays questions in GUI)
%
% function z = askparam(ask, prompt, zLow, zUpp, zDef, UP, Itext,UPpos)
%
% If ask > 0 asks question
% If ask == 0 uses default value zDef
% If ask < 0 uses UP(1) is nonempty, otherwise zDef
%
% INPUT:
% prompt    Prompt string.             Default: Give real value
% zLow      Lower bound on real value. Default -Inf
% zUpp      Upper bound on real value. Default  Inf
% zDef      Default value for z if defined, otherwise not used.
% ask       If ask > 10 GUI question, otherwise in MATLAB window
% UP        User parameters, UP(1) is used as default value (if ask < 0 )
% Itext     Extra instruction text, displayed before question.
% UPpos     Position in UP vector, default 1
%
% OUTPUT:
% z         Real value in [zLow,zUpp], input from user choice

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Feb 19, 1998.   Last modified March 23, 2003.

function z = askparam(ask, prompt, zLow, zUpp, zDef, UP, Itext,UPpos)

if nargin < 8
   UPpos=[];
   if nargin < 7
      Itext=[];
      if nargin < 6
         UP=[];
         if nargin < 5
            zDef=[];
            if nargin < 4
               zUpp=[];
               if nargin < 3
                  zLow=[];
                  if nargin < 2
                     prompt=[];
                  end
               end
            end
         end
      end
   end
end
if isempty(UPpos),  UPpos=1;                   end
if isempty(zUpp),   zUpp=Inf;                  end
if isempty(zLow),   zLow=-Inf;                 end
if isempty(prompt), prompt='Give real value '; end

if ask>=1
   if length(UP) >= UPpos 
      z=UP(UPpos);
   else
      z=zDef;
   end
   %if isempty(z), disp('IS EMPTY'); end
   T=str2mat(Itext,' ',['Current value = ' num2str(z)]);
   z=inputR(prompt,zLow,zUpp,z,ask,T); 
elseif ask==0
   z = zDef;
else
   if length(UP) >= UPpos 
      z=UP(UPpos);
   else
      z=zDef;
   end
end