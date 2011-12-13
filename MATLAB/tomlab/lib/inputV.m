% InputV prompts the user for a vector with length n 
% with values in the interval [Rlow,Rupp]
% If dupl is set, allow 1 element answer, which is duplicated
%
% If the user gives the wrong answer, an error message, the interval, and
% the wanted vector length is displayed.
%
% function v = inputV(prompt, n, Rlow, Rupp, dupl, ask, Itext)
%
% INPUT:
% prompt    Prompt string.             Default: Give vector
% n         Length of vector.          Default    2
% Rlow      Lower bound on real value. Default -Inf
% Rupp      Upper bound on real value. Default  Inf
% dupl      If dupl==1, allow 1 element which is duplicated. Default dupl=0.
% ask       If ask > 10 GUI question, otherwise in MATLAB window
% Itext     Extra instruction text, displayed before question.
%
% OUTPUT:
% v         Vector with values in [Rlow,Rupp]. v is ALWAYS a column vector.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1997-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.9.0$
% Written Mar 26, 1997.  Last modified Aug 01, 2005.

function v = inputV(prompt, n, Rlow, Rupp, dupl, ask, Itext)

if nargin < 7
   Itext=[];
   if nargin < 6
      ask=[];
      if nargin < 5
         dupl=[];
         if nargin < 4
            Rupp=[];
            if nargin < 3
               Rlow=[];
               if nargin < 2
                  n=[];
                  if nargin < 1
                     prompt=[];
                  end
               end
            end
         end
      end
   end
end

global question answer instruction

if isempty(ask),    ask=1;                 end
if isempty(Rupp),   Rupp=Inf;              end
if isempty(Rlow),   Rlow=-Inf;             end
if isempty(prompt), prompt='Give vector '; end
if isempty(n),      n=2;                   end
if isempty(dupl),   dupl=0;                end

if ask > 10
   question=prompt;
end
instruction=Itext;

while 1
   if ask > 10
      %feval('tomGUI','askquest');
      askgui('askquest');
      answer=deblank(answer);
      if isempty(answer)
         v=[];
      else
         v=str2num(answer);
      end
   else
      disp(instruction);
      fprintf('\n');
      v=input(prompt,'s');
      v=str2num(v);
   end
   if dupl & length(v) == 1 
      if v >= Rlow & v <= Rupp
         v=v*ones(n,1);
         return;
      else
         instruction=['ERROR! The value is outside the interval [' ...
              num2str(Rlow) ',' num2str(Rupp) '].'];
      end
   elseif length(v) == n
      % Check limits
      v=v(:);
      if all(v >= Rlow) & all(v <= Rupp), return; end

      instruction=...
        ['ERROR! The vector has one or more values outside the interval [' ...
         num2str(Rlow) ',' num2str(Rupp) '].'];
   else
      instruction=['Your answer has not the length ' num2str(n) '.'];
      if isempty(v)
         instruction=str2mat(instruction,'Your answer is empty!');
      elseif ischar(v)
         instruction=str2mat(instruction,...
           'Your answer is a string, not numeric!');
      end
   end
end

% MODIFICATION LOG:
%
% 981208 hkh  Safequard, test empty on GUI answer, otherwise crash if user 
%             gives ENTER and then OK. The anwer may have some blanks in it.
% 990618 hkh  Change to use tomlabGUI
% 990726 hkh  Use askgui instead of tomlabGUI
% 990911 hkh  Use str2num instead of eval
% 000927 hkh  Serious bug, must use input to get the string when not GUI
% 050801 med  isstr changed to ischar