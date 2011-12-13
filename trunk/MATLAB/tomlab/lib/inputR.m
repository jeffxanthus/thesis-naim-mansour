% InputR prompts the user for one real value z in the interval [zLow,zUpp]
%
% If the user gives the wrong answer, the interval is displayed.
% If the user gives a vector, 1st element is taken as the response
% If zDef is defined, zDef is used as default value if user gives an
% empty answer. No check is made on the zDef value.
%
% function z = inputR(prompt, zLow, zUpp, zDef, ask, Itext)
%
% INPUT:
% prompt    Prompt string.             Default: Give real value
% zLow      Lower bound on real value. Default -Inf
% zUpp      Upper bound on real value. Default  Inf
% zDef      Default value for z if defined, otherwise not used.
% ask       If ask > 10 GUI question, otherwise in MATLAB window
% Itext     Extra instruction text, displayed before question.
%
% OUTPUT:
% z         Real value in [zLow,zUpp], input from user choice
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1997-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Mar 25, 1997.  Last modified Sep 10, 1999.

function z = inputR(prompt, zLow, zUpp, zDef, ask, Itext)

if nargin < 6
   Itext=[];
   if nargin < 5
      ask=[];
      if nargin < 4
         zDef=[];
         if nargin < 3
            zUpp=[];
            if nargin < 2
               zLow=[];
               if nargin < 1
                  prompt=[];
               end
            end
         end
      end
   end
end

global question answer instruction

if isempty(ask),    ask=1;                     end
if isempty(zUpp),   zUpp=Inf;                  end
if isempty(zLow),   zLow=-Inf;                 end
if isempty(prompt), prompt='Give real value '; end

if ask > 10
   question=prompt;
end
instruction=Itext;

if ~isempty(zDef)
   zDef=min(max(zDef,zLow),zUpp);
end

while 1
   if ask > 10
      %feval('tomGUI','askquest');
      askgui('askquest');
      answer=deblank(answer);
      if isempty(answer)
         z=[];
      else
         z=str2num(answer);
      end
   else
      disp(instruction);
      fprintf('\n');
      z=input(prompt);
      % Using inputdlg gives no error checking
      % Also cannot handle several instruction rows
      % if size(instruction,1) > 1
      %    disp(instruction);
      % end
      % zz=inputdlg(prompt,instruction(size(instruction,1),:),1);
      % z=char(zz(1));
      % z=deblank(z);
      % if ~isempty(z)
      %    z=str2num(z);
      % else
      %    z=[];
      % end
   end
   if isempty(z) & ~isempty(zDef)
      z=zDef;
      return
   elseif isempty(z) 
      instruction=['Your empty answer is not in the interval ' ...
           num2str(zLow) ', ' num2str(zUpp)];
   else
      z=z(1);
      if z >= zLow & z <= zUpp, return; end

      instruction=['Your answer ' num2str(z) ' is not in the interval ' ...
           num2str(zLow) ', ' num2str(zUpp)];
   end
end

% MODIFICATION LOG
%
% 981208 hkh  Safequard, test empty on GUI answer, otherwise crash if user
%             gives ENTER and then OK.
% 981209 hkh  Put zDef in [zLow,zUpp]
% 990618 hkh  Change to use tomlabGUI
% 990726 hkh  Use askgui instead of tomlabGUI
% 990911 hkh  Use str2num instead of eval