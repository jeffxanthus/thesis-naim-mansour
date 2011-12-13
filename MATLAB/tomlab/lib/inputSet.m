% InputSet prompts the user for 1 integer value among the set SetVals.
%
% If the user gives the wrong answer, the set is displayed.
% If the user gives a vector, 1st element is taken as the response
%
% function k = inputSet(prompt, SetVals, ask, Itext)
%
% INPUT:
% prompt    Prompt string
% SetVals   Set of integers which user has to choose from
% ask       If ask > 10 GUI question, otherwise in MATLAB window
% Itext     Extra instruction text, displayed before question.
%
% OUTPUT:
% k         User choice

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1997-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Mar 25, 1997.  Last modified Sep 11, 1999.

function k = inputSet(prompt, SetVals, ask, Itext)

if nargin < 4
   Itext=[];
   if nargin < 3
      ask=[];
   end
end

global question answer instruction

if isempty(ask), ask=1; end

if ask > 10
   question=prompt;
   set=['[' num2str(SetVals(1))];
   for i=2:length(SetVals)
      set = [set ' ' num2str(SetVals(i))];
   end
   set = [set ']'];
end

instruction=Itext;

while 1
   if ask > 10
      %feval('tomGUI','askquest');
      askgui('askquest');
      answer=deblank(answer);
      if isempty(answer)
         k=[];
      else
         k=str2num(answer);
      end
   else
      disp(instruction);
      fprintf('\n');
      k=input(prompt);
   end
   if ~isempty(k)
      k=k(1);
      if any(k==SetVals), return, end

      instruction=str2mat(...
          ['Your answer ' num2str(k) ' is wrong. Choose between:'],set);
   else
      instruction=str2mat('Your empty answer is wrong. Choose between:',set);
   end
end

% MODIFICATION LOG:
%
% 981208 hkh  Safequard, test empty on GUI answer, otherwise crash if user 
%             gives ENTER and then OK. May also be some blanks in the answer.
% 990618 hkh  Change to use tomlabGUI
% 990726 hkh  Use askgui instead of tomlabGUI
% 990911 hkh  Use str2num instead of eval